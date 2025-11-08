package sunyu.util;

import java.util.ArrayList;
import java.util.List;

import org.geotools.geometry.jts.JTS;
import org.geotools.referencing.CRS;
import org.geotools.referencing.crs.DefaultGeographicCRS;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.GeometryFactory;
import org.opengis.referencing.crs.CoordinateReferenceSystem;
import org.opengis.referencing.operation.MathTransform;

import cn.hutool.log.Log;
import cn.hutool.log.LogFactory;
import sunyu.util.pojo.SplitRoadResult;
import sunyu.util.pojo.TrackPoint;

/**
 * GIS工具类，封装轨迹轮廓构建、道路分段、坐标转换、拓扑判断与面积计算。
 * 
 * 
 * @author SunYu
 */
public class GisUtil implements AutoCloseable {
    private final Log log = LogFactory.get();
    // 配置参数，包含各种常量和默认值
    private final Config config;

    /**
     * 创建构建器实例（Builder）。
     * 通过构建器配置参数后再构建 `GisUtil`，避免直接实例化带来的不完整/不一致配置。
     *
     * @return 用于配置并构建 `GisUtil` 的构建器
     */
    public static Builder builder() {
        return new Builder();
    }

    /**
     * 私有构造函数，通过 Builder 创建实例。
     * 保持不可直接实例化，集中初始化日志与配置引用。
     *
     * @param config 内部配置对象，包含常量、默认值与缓存实例
     */
    private GisUtil(Config config) {
        log.info("[构建{}] 开始", this.getClass().getSimpleName());
        log.info("[构建{}] 结束", this.getClass().getSimpleName());
        this.config = config;
    }

    /**
     * 内部配置类。
     * <p>
     * 持有常量、默认参数、线程安全缓存与几何工厂；供 {@link GisUtil} 使用。
     * 该类包含所有GIS处理过程中需要的配置参数，确保线程安全的参数访问。
     * </p>
     */
    private static class Config {
        // 几何工厂
        public final GeometryFactory geometryFactory = new GeometryFactory();
        // 空几何集合 - 使用空数组参数更明确和安全
        public final Geometry EMPTYGEOM = geometryFactory.createGeometryCollection(new Geometry[0]);
        // 空WKT字符串
        public final String EMPTY_WKT = EMPTYGEOM.toText();

        // 最大速度（单位：米/秒），我们认为农机在田间作业，最大速度不会超过20米/秒
        public final double MAX_SPEED = 20.0;
    }

    /**
     * 构建器类，用于构建GisUtil实例。
     * <p>
     * 实现构建器模式，允许灵活配置和创建GisUtil对象，提供流式API进行参数设置。
     * </p>
     */
    public static class Builder {
        // 配置对象，包含各种常量和默认值
        private Config config = new Config();

        /**
         * 构建GisUtil实例。
         * <p>
         * 使用预配置的参数创建GisUtil对象，完成所有必要的初始化工作。
         * </p>
         * 
         * @return 已初始化完成的GisUtil对象，可直接使用
         */
        public GisUtil build() {
            return new GisUtil(config);
        }

    }

    /**
     * 关闭资源（AutoCloseable）
     * 重置转换对象以释放资源。
     */
    @Override
    public void close() {
        log.info("[销毁{}] 开始", this.getClass().getSimpleName());
        log.info("[销毁{}] 结束", this.getClass().getSimpleName());
    }

    /**
     * 将WGS84坐标点转换为高斯投影坐标点
     * 支持全球范围，根据经度自动选择合适的高斯投影带
     * 
     * @param point WGS84坐标系下的轨迹点
     * @return 高斯投影坐标系下的轨迹点
     */
    private TrackPoint wgs84PointTransformToGaussPoint(TrackPoint point) {
        try {
            // 获取经度和纬度
            double longitude = point.getLon();
            double latitude = point.getLat();

            log.debug("开始坐标转换：原始经度={}, 纬度={}", longitude, latitude);

            // 计算高斯投影带号 (6度分带)
            // 全球范围: 经度-180到180，对应带号1-60
            int zone = (int) Math.floor((longitude + 180) / 6) + 1;
            log.debug("计算投影带号：zone={}", zone);

            // 计算中央经线
            double centralMeridian = (zone - 1) * 6 - 180 + 3;
            log.debug("计算中央经线：centralMeridian={}", centralMeridian);

            // 验证投影带号的合理性
            if (zone < 1 || zone > 60) {
                log.warn("投影带号超出合理范围：zone={}，经度={}", zone, longitude);
            }

            // 创建WGS84坐标参考系统
            CoordinateReferenceSystem wgs84CRS = DefaultGeographicCRS.WGS84;

            // 创建高斯投影坐标参考系统
            // 全球支持模式：假东距包含带号信息，便于识别投影带
            // 格式：zone × 1000000 + 500000（如49带 = 49500000米）
            double falseEasting = zone * 1000000.0 + 500000.0;
            log.debug("计算假东距：falseEasting={}", falseEasting);
            String wkt = String.format(
                    "PROJCS[\"WGS_1984_Gauss_Kruger_Zone_%d\", " +
                            "GEOGCS[\"GCS_WGS_1984\", " +
                            "DATUM[\"D_WGS_1984\", " +
                            "SPHEROID[\"WGS_1984\",6378137.0,298.257223563]], " +
                            "PRIMEM[\"Greenwich\",0.0], " +
                            "UNIT[\"Degree\",0.0174532925199433]], " +
                            "PROJECTION[\"Transverse_Mercator\"], " +
                            "PARAMETER[\"False_Easting\",%.1f], " +
                            "PARAMETER[\"False_Northing\",0.0], " +
                            "PARAMETER[\"Central_Meridian\",%.1f], " +
                            "PARAMETER[\"Scale_Factor\",1.0], " +
                            "PARAMETER[\"Latitude_Of_Origin\",0.0], " +
                            "UNIT[\"Meter\",1.0]]",
                    zone, falseEasting, centralMeridian);

            log.debug("WKT定义：{}", wkt);

            CoordinateReferenceSystem gaussCRS = CRS.parseWKT(wkt);

            // 创建坐标转换
            MathTransform transform = CRS.findMathTransform(wgs84CRS, gaussCRS, true);

            // 执行坐标转换
            Coordinate sourceCoord = new Coordinate(longitude, latitude);
            Coordinate targetCoord = new Coordinate();

            JTS.transform(sourceCoord, targetCoord, transform);

            log.debug("转换结果：X={}, Y={}", targetCoord.x, targetCoord.y);

            // 验证转换结果的合理性（全球范围）
            // X坐标：最小1带=150万米，最大60带=6050万米，加上实际坐标范围±350万米
            // 合理范围：100万-6400万米
            // Y坐标：全球纬度范围对应约±1000万米
            if (targetCoord.x < 1000000 || targetCoord.x > 64000000 || targetCoord.y < -10000000
                    || targetCoord.y > 10000000) {
                log.warn("转换结果超出全球合理范围：X={}, Y={}, zone={}", targetCoord.x, targetCoord.y, zone);
            }

            // 创建新的轨迹点，保持原有属性，只更新坐标
            TrackPoint result = new TrackPoint();
            result.setTime(point.getTime());
            result.setLon(targetCoord.x); // X坐标（东向）
            result.setLat(targetCoord.y); // Y坐标（北向）
            result.setSpeed(point.getSpeed());
            result.setDirection(point.getDirection());

            log.debug("坐标转换成功完成");
            return result;

        } catch (Exception e) {
            log.error("坐标转换失败：经度={}, 纬度={}, 错误={}", point.getLon(), point.getLat(), e.getMessage());
            return null;
        }
    }

    public SplitRoadResult splitRoad(List<TrackPoint> wgs84Points, double totalWidthM) {
        SplitRoadResult result = new SplitRoadResult();
        result.setOutline(config.EMPTYGEOM);
        result.setWkt(config.EMPTY_WKT);

        log.debug("参数信息：轨迹点数量={} 总宽度={}米", wgs84Points.size(), totalWidthM);
        log.debug("过滤时间为空、经纬度异常、经纬度为0、速度为0、方向角异常的轨迹点");
        List<TrackPoint> filteredWgs84Points = new ArrayList<>();
        for (TrackPoint point : wgs84Points) {
            if (point.getTime() != null// 定位时间不能为空
                    && point.getLon() >= -180 && point.getLon() <= 180// 经度必须在[-180,180]之间
                    && point.getLat() >= -90 && point.getLat() <= 90// 纬度必须在[-90,90]之间
                    && point.getLon() != 0 && point.getLat() != 0// 经度和纬度不能为0
                    && point.getSpeed() > 0// 速度必须大于0
                    && point.getDirection() >= 0 && point.getDirection() <= 360) {// 方向角必须在[0,360]之间
                filteredWgs84Points.add(point);
            } else {
                log.warn("被过滤点位信息，时间：{} 经度：{} 纬度：{} 速度：{} 方向角：{}", point.getTime(), point.getLon(), point.getLat(),
                        point.getSpeed(), point.getDirection());
            }
        }
        log.debug("过滤后的轨迹点数量：{}", filteredWgs84Points.size());

        log.debug("按定位时间升序排序");
        filteredWgs84Points.sort((p1, p2) -> p1.getTime().compareTo(p2.getTime()));

        log.debug("按照速度切分成多段轨迹，只要速度超过{}米/秒，就拆分轨迹段", config.MAX_SPEED);
        List<List<TrackPoint>> wgs84PointsSegments = new ArrayList<>();
        List<TrackPoint> currentSegment = new ArrayList<>();
        for (TrackPoint point : filteredWgs84Points) {
            if (point.getSpeed() >= config.MAX_SPEED) {
                if (!currentSegment.isEmpty()) {
                    wgs84PointsSegments.add(currentSegment);
                    currentSegment = new ArrayList<>();
                }
            }
            currentSegment.add(point);
        }
        if (!currentSegment.isEmpty()) {
            wgs84PointsSegments.add(currentSegment);
        }
        log.debug("速度切分后的轨迹段数量：{}", wgs84PointsSegments.size());

        log.debug("过滤掉段内小于6个点的轨迹段");
        wgs84PointsSegments.removeIf(segment -> segment.size() < 6);
        log.debug("过滤后的轨迹段数量：{}", wgs84PointsSegments.size());

        for (List<TrackPoint> wgs84PointsSegment : wgs84PointsSegments) {
            log.debug("将轨迹段中的WGS84坐标转换为高斯投影坐标 轨迹段点数量：{}", wgs84PointsSegment.size());
            List<TrackPoint> gaussPoints = new ArrayList<>();
            for (TrackPoint point : wgs84PointsSegment) {
                TrackPoint gaussPoint = wgs84PointTransformToGaussPoint(point);
                if (gaussPoint != null) {
                    gaussPoints.add(gaussPoint);
                }
            }
            log.debug("转换后的轨迹段点数量：{}", gaussPoints.size());
        }

        return result;
    }

}