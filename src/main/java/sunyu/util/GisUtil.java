package sunyu.util;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.ConcurrentHashMap;

import org.geotools.geometry.jts.JTS;
import org.geotools.referencing.CRS;
import org.geotools.referencing.crs.DefaultGeographicCRS;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Envelope;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.LineString;
import org.opengis.referencing.crs.CoordinateReferenceSystem;
import org.opengis.referencing.operation.MathTransform;

import cn.hutool.log.Log;
import cn.hutool.log.LogFactory;
import sunyu.util.pojo.OutlinePart;
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

        // WGS84坐标参考系统（全局常量，避免重复创建）
        public final CoordinateReferenceSystem WGS84_CRS = DefaultGeographicCRS.WGS84;

        // 高斯投影CRS缓存（线程安全），key格式："zone_falseEasting_centralMeridian"
        public final ConcurrentHashMap<String, CoordinateReferenceSystem> gaussCRSCache = new ConcurrentHashMap<>();

        // WGS84到高斯投影的坐标转换缓存（线程安全），key格式："zone_falseEasting_centralMeridian"
        public final ConcurrentHashMap<String, MathTransform> wgs84ToGaussTransformCache = new ConcurrentHashMap<>();

        // 高斯投影到WGS84的坐标转换缓存（线程安全），key格式："zone_falseEasting_centralMeridian"
        public final ConcurrentHashMap<String, MathTransform> gaussToWgs84TransformCache = new ConcurrentHashMap<>();

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
     * 从缓存获取或创建高斯投影CRS
     * 
     * @param zone            投影带号
     * @param falseEasting    假东距
     * @param centralMeridian 中央经线
     * @return 高斯投影CRS对象，如果创建失败返回null
     */
    private CoordinateReferenceSystem getOrCreateGaussCRS(int zone, double falseEasting, double centralMeridian) {
        // 构建缓存key：zone_falseEasting_centralMeridian
        String cacheKey = String.format("%d_%.1f_%.1f", zone, falseEasting, centralMeridian);

        // 从缓存获取或创建高斯投影CRS
        return config.gaussCRSCache.computeIfAbsent(cacheKey, key -> {
            try {
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

                log.debug("高斯投影CRS WKT定义：{}", wkt);
                return CRS.parseWKT(wkt);
            } catch (Exception e) {
                log.warn("解析WKT失败：zone={}, falseEasting={}, centralMeridian={}, 错误={}",
                        zone, falseEasting, centralMeridian, e.getMessage());
                return null;
            }
        });
    }

    /**
     * 将WGS84坐标点转换为高斯投影坐标点
     * 支持全球范围，根据经度自动选择合适的高斯投影带
     * 
     * @param wgs84Point WGS84坐标系下的轨迹点
     * @return 高斯投影坐标系下的轨迹点
     */
    private TrackPoint wgs84PointTransformToGaussPoint(TrackPoint wgs84Point) {
        try {
            // 获取经度和纬度
            double longitude = wgs84Point.getLon();
            double latitude = wgs84Point.getLat();

            log.trace("=== 开始WGS84到高斯投影转换 ===");
            log.trace("原始坐标：经度={}, 纬度={}", longitude, latitude);

            // 计算高斯投影带号 (6度分带)
            // 全球范围: 经度-180到180，对应带号1-60
            int zone = (int) Math.floor((longitude + 180) / 6) + 1;
            log.trace("计算投影带号：zone={}", zone);

            // 计算中央经线
            double centralMeridian = (zone - 1) * 6 - 180 + 3;
            log.trace("计算中央经线：centralMeridian={}", centralMeridian);

            // 验证投影带号的合理性
            if (zone < 1 || zone > 60) {
                log.warn("投影带号超出合理范围：zone={}，经度={}", zone, longitude);
                return null;
            }

            // 创建高斯投影坐标参考系统（使用缓存避免重复解析WKT）
            // 全球支持模式：假东距包含带号信息，便于识别投影带
            // 格式：zone × 1000000 + 500000（如49带 = 49500000米）
            double falseEasting = zone * 1000000.0 + 500000.0;
            log.trace("计算假东距：falseEasting={}", falseEasting);

            // 从缓存获取或创建高斯投影CRS
            CoordinateReferenceSystem gaussCRS = getOrCreateGaussCRS(zone, falseEasting, centralMeridian);

            // 构建缓存key（用于transform缓存）
            String cacheKey = String.format("%d_%.1f_%.1f", zone, falseEasting, centralMeridian);

            // 从缓存获取或创建WGS84到高斯投影的坐标转换
            MathTransform transform = config.wgs84ToGaussTransformCache.computeIfAbsent(cacheKey, key -> {
                try {
                    return CRS.findMathTransform(config.WGS84_CRS, gaussCRS, true);
                } catch (Exception e) {
                    log.warn("创建坐标转换失败：zone={}, falseEasting={}, centralMeridian={}, 错误={}",
                            zone, falseEasting, centralMeridian, e.getMessage());
                    return null;
                }
            });

            // 执行坐标转换
            Coordinate sourceCoord = new Coordinate(longitude, latitude);
            Coordinate targetCoord = new Coordinate();

            JTS.transform(sourceCoord, targetCoord, transform);

            log.trace("转换结果：X={}, Y={}", targetCoord.x, targetCoord.y);

            // 验证转换结果的合理性（全球范围）
            // X坐标：最小1带=150万米，最大60带=6050万米，加上实际坐标范围±350万米
            // 合理范围：50万-6400万米（更宽松的范围，允许边界情况）
            // Y坐标：全球纬度范围对应约±1000万米
            if (targetCoord.x < 500000 || targetCoord.x > 64000000 || targetCoord.y < -10000000
                    || targetCoord.y > 10000000) {
                log.warn("转换结果超出全球合理范围：X={}, Y={}, zone={}", targetCoord.x, targetCoord.y, zone);
                return null;
            }

            // 创建新的轨迹点，保持原有属性，只更新坐标
            TrackPoint result = new TrackPoint();
            result.setTime(wgs84Point.getTime());
            result.setLon(targetCoord.x); // X坐标（东向）
            result.setLat(targetCoord.y); // Y坐标（北向）
            result.setSpeed(wgs84Point.getSpeed());
            result.setDirection(wgs84Point.getDirection());

            log.trace("=== WGS84到高斯投影转换完成 ===");
            return result;
        } catch (Exception e) {
            log.warn("坐标转换失败：经度={}, 纬度={}, 错误={}", wgs84Point.getLon(), wgs84Point.getLat(), e.getMessage());
            return null;
        }
    }

    private Geometry gaussGeometryToWgs84Geometry(Geometry gaussGeometry) {
        try {
            // 获取几何图形的边界信息来反推高斯投影参数
            Envelope env = gaussGeometry.getEnvelopeInternal();
            double centerX = (env.getMinX() + env.getMaxX()) / 2.0;

            log.debug("高斯投影几何到WGS84投影几何转换开始");
            log.trace("高斯几何边界：MinX={}, MaxX={}, MinY={}, MaxY={}",
                    env.getMinX(), env.getMaxX(), env.getMinY(), env.getMaxY());
            log.trace("几何中心X坐标：centerX={}", centerX);

            // 验证高斯投影坐标的合理性
            // 高斯投影X坐标合理范围：50万-6400万米（对应全球1-60带，更宽松的范围）
            // 高斯投影Y坐标合理范围：±1000万米（对应全球纬度范围）
            if (env.getMinX() < 500000 || env.getMaxX() > 64000000 ||
                    env.getMinY() < -10000000 || env.getMaxY() > 10000000) {
                log.warn("高斯投影坐标超出合理范围：MinX={}, MaxX={}, MinY={}, MaxY={}",
                        env.getMinX(), env.getMaxX(), env.getMinY(), env.getMaxY());
                return config.EMPTYGEOM;
            }

            // 智能确定投影带号策略
            // 1. 首先尝试根据几何中心反推带号
            int zone = (int) Math.floor(centerX / 1000000.0);
            log.trace("反推计算的投影带号：zone={}, centerX={}", zone, centerX);

            // 2. 如果几何范围跨越多个投影带（宽度超过100万米），使用更保守的策略
            double geometryWidth = env.getMaxX() - env.getMinX();
            if (geometryWidth > 1000000) {
                log.warn("几何图形宽度{}米，可能跨越多个投影带，使用保守策略", geometryWidth);
                // 对于宽几何，优先使用中间带号，并验证整个几何范围
                if (zone < 1 || zone > 60) {
                    log.warn("无法确定合适的投影带号：centerX={}", centerX);
                    return config.EMPTYGEOM;
                }
            } else if (zone < 1 || zone > 60) {
                log.warn("反推的投影带号不合理：zone={}, centerX={}", zone, centerX);
                // 尝试备用策略：如果zone不合理，可能是假东距计算问题，尝试重新计算
                log.trace("尝试备用策略重新计算投影带号");
                int backupZone = (int) Math.floor((centerX - 500000.0) / 1000000.0);
                log.trace("备用策略计算的投影带号：zone={}", backupZone);
                if (backupZone < 1 || backupZone > 60) {
                    log.warn("备用策略仍然无法确定合适的投影带号：centerX={}", centerX);
                    return config.EMPTYGEOM;
                }
                zone = backupZone;
            }

            // 计算中央经线（与wgs84PointTransformToGaussPoint方法保持一致）
            double centralMeridian = (zone - 1) * 6 - 180 + 3;
            double falseEasting = zone * 1000000.0 + 500000.0;
            log.trace("计算投影参数：zone={}, centralMeridian={}, falseEasting={}", zone, centralMeridian, falseEasting);

            // 从缓存获取或创建高斯投影CRS
            CoordinateReferenceSystem gaussCRS = getOrCreateGaussCRS(zone, falseEasting, centralMeridian);

            if (gaussCRS == null) {
                log.warn("无法获取高斯投影CRS：zone={}", zone);
                return config.EMPTYGEOM;
            }

            // 构建缓存key（用于transform缓存）
            String cacheKey = String.format("%d_%.1f_%.1f", zone, falseEasting, centralMeridian);

            // 从缓存获取或创建高斯投影到WGS84的坐标转换
            final int finalZone = zone;
            MathTransform transform = config.gaussToWgs84TransformCache.computeIfAbsent(cacheKey, key -> {
                try {
                    return CRS.findMathTransform(gaussCRS, config.WGS84_CRS, true);
                } catch (Exception e) {
                    log.warn("创建坐标转换失败：zone={}, falseEasting={}, centralMeridian={}, 错误={}",
                            finalZone, falseEasting, centralMeridian, e.getMessage());
                    return null;
                }
            });

            if (transform == null) {
                log.warn("无法获取坐标转换：zone={}", finalZone);
                return config.EMPTYGEOM;
            }

            // 执行坐标转换（逆向：高斯投影 -> WGS84）
            log.trace("zone={}, 几何点数={}", finalZone, gaussGeometry.getNumPoints());
            Geometry wgs84Geometry = JTS.transform(gaussGeometry, transform);

            // 验证转换后的WGS84坐标合理性
            Envelope wgs84Env = wgs84Geometry.getEnvelopeInternal();
            if (wgs84Env.getMinX() < -180 || wgs84Env.getMaxX() > 180 ||
                    wgs84Env.getMinY() < -90 || wgs84Env.getMaxY() > 90) {
                log.warn("转换后的WGS84坐标超出合理范围：MinLon={}, MaxLon={}, MinLat={}, MaxLat={}",
                        wgs84Env.getMinX(), wgs84Env.getMaxX(), wgs84Env.getMinY(), wgs84Env.getMaxY());
                return config.EMPTYGEOM;
            }
            log.debug("高斯投影几何到WGS84投影几何转换完成");
            return wgs84Geometry;
        } catch (Exception e) {
            log.warn("高斯投影几何到WGS84投影几何转换失败：错误={}", e.getMessage());
            return config.EMPTYGEOM;
        }
    }

    /**
     * 将高斯投影的几何图形转换为WGS84坐标系的WKT字符串
     * 
     * @param gaussGeometry 高斯投影坐标系下的几何图形
     * @return WGS84坐标系下的WKT字符串
     */
    private String gaussGeometryToWgs84WKT(Geometry gaussGeometry) {
        log.debug("高斯投影几何转换为WGS84投影WKT字符串开始");
        Geometry wgs84Geometry = gaussGeometryToWgs84Geometry(gaussGeometry);
        log.debug("高斯投影几何转换为WGS84投影WKT字符串完成");
        return wgs84Geometry.toText();
    }

    public SplitRoadResult splitRoad(List<TrackPoint> wgs84Points, double totalWidthM) {
        SplitRoadResult result = new SplitRoadResult();
        result.setOutline(config.EMPTYGEOM);
        result.setWkt(config.EMPTYGEOM.toText());

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
                log.trace("被过滤点位信息，时间：{} 经度：{} 纬度：{} 速度：{} 方向角：{}", point.getTime(), point.getLon(), point.getLat(),
                        point.getSpeed(), point.getDirection());
            }
        }
        log.debug("过滤后的轨迹点数量：{}", filteredWgs84Points.size());
        if (filteredWgs84Points.size() < 6) {
            log.warn("过滤后的轨迹点数量小于6个，无法进行拆分");
            return result;
        }

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

        List<OutlinePart> outlineParts = new ArrayList<>();
        for (List<TrackPoint> wgs84PointsSegment : wgs84PointsSegments) {
            if (wgs84PointsSegment.size() < 6) {
                log.warn("轨迹段点数量小于6个，无法进行拆分");
                continue;
            }

            log.debug("将轨迹段中的WGS84坐标转换为高斯投影坐标 轨迹段点数量：{}", wgs84PointsSegment.size());
            List<TrackPoint> gaussPoints = new ArrayList<>();
            for (TrackPoint point : wgs84PointsSegment) {
                TrackPoint gaussPoint = wgs84PointTransformToGaussPoint(point);
                if (gaussPoint != null) {
                    gaussPoints.add(gaussPoint);
                }
            }
            log.debug("转换后的轨迹段点数量：{}", gaussPoints.size());
            if (gaussPoints.size() < 6) {
                log.warn("转换后的轨迹段点数量小于6个，无法进行拆分");
                continue;
            }

            log.debug("使用高斯投影的轨迹点进行线缓冲，左右缓冲总宽度：{}米", totalWidthM);
            try {
                // 1. 创建线几何：将高斯投影点转换为Coordinate数组
                Coordinate[] gaussCoordinates = new Coordinate[gaussPoints.size()];
                for (int i = 0; i < gaussPoints.size(); i++) {
                    TrackPoint point = gaussPoints.get(i);
                    gaussCoordinates[i] = new Coordinate(point.getLon(), point.getLat());
                }

                // 2. 创建LineString
                LineString lineString = config.geometryFactory.createLineString(gaussCoordinates);
                log.debug("创建线几何成功，点数：{}", gaussCoordinates.length);

                // 3. 执行缓冲操作：总宽度的一半作为缓冲距离
                double bufferDistance = totalWidthM / 2.0;
                Geometry bufferedGeometry = lineString.buffer(bufferDistance);
                log.debug("线缓冲成功，缓冲距离：{}米，结果几何类型：{}", bufferDistance, bufferedGeometry.getGeometryType());

                OutlinePart outlinePart = new OutlinePart();
                outlinePart.setOutline(bufferedGeometry);
                outlinePart.setWkt(gaussGeometryToWgs84WKT(bufferedGeometry));
                outlinePart.setTrackPoints(wgs84PointsSegment);
                outlinePart.setStartTime(wgs84PointsSegment.get(0).getTime());
                outlinePart.setEndTime(wgs84PointsSegment.get(wgs84PointsSegment.size() - 1).getTime());
                outlineParts.add(outlinePart);
            } catch (Exception e) {
                log.warn("线缓冲失败：总宽度={}米，错误={}", totalWidthM, e.getMessage());
            }
        }

        if (outlineParts.isEmpty()) {
            log.warn("没有成功创建任何区块，无法进行合并");
            return result;
        }

        log.debug("合并所有高斯投影缓冲几何");
        Geometry unionGeometry = config.geometryFactory
                .createGeometryCollection(outlineParts.stream().map(OutlinePart::getOutline).toArray(Geometry[]::new))
                .union();
        log.debug("合并后的几何类型：{}", unionGeometry.getGeometryType());

        // 4. 设置结果
        result.setOutline(unionGeometry);
        result.setWkt(gaussGeometryToWgs84WKT(unionGeometry));
        result.setParts(outlineParts);

        return result;
    }

}