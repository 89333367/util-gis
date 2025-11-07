package sunyu.util;

import java.time.Duration;
import java.time.LocalDateTime;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.concurrent.ConcurrentHashMap;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import org.geotools.geometry.DirectPosition2D;
import org.geotools.geometry.jts.JTS;
import org.geotools.referencing.CRS;
import org.geotools.referencing.crs.DefaultGeographicCRS;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Envelope;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.LineString;
import org.locationtech.jts.geom.MultiPolygon;
import org.locationtech.jts.geom.Point;
import org.locationtech.jts.geom.Polygon;
import org.locationtech.jts.geom.TopologyException;
import org.locationtech.jts.geom.prep.PreparedGeometry;
import org.locationtech.jts.geom.prep.PreparedGeometryFactory;
import org.locationtech.jts.io.ParseException;
import org.locationtech.jts.io.WKTReader;
import org.locationtech.jts.operation.buffer.BufferOp;
import org.locationtech.jts.operation.buffer.BufferParameters;
import org.locationtech.jts.operation.overlay.OverlayOp;
import org.locationtech.jts.operation.overlay.snap.SnapIfNeededOverlayOp;
import org.locationtech.jts.operation.union.UnaryUnionOp;
import org.locationtech.jts.simplify.DouglasPeuckerSimplifier;
import org.opengis.referencing.crs.CoordinateReferenceSystem;
import org.opengis.referencing.operation.MathTransform;

import cn.hutool.core.collection.CollUtil;
import cn.hutool.core.util.StrUtil;
import cn.hutool.log.Log;
import cn.hutool.log.LogFactory;
import sunyu.util.pojo.CoordinatePoint;
import sunyu.util.pojo.OutlinePart;
import sunyu.util.pojo.SplitRoadResult;
import sunyu.util.pojo.TrackPoint;
import sunyu.util.pojo.WktIntersectionResult;

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
        // 记录工具类构建开始日志
        log.info("[构建{}] 开始", this.getClass().getSimpleName());
        // 其他初始化语句（预留扩展点）
        // 记录工具类构建结束日志
        log.info("[构建{}] 结束", this.getClass().getSimpleName());
        // 保存配置参数引用
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
        // WGS84坐标系的EPSG代码，用于定义地理坐标系统
        private final String WGS84 = "EPSG:4326";

        // 地球半径（米），用于Haversine公式计算两点间距离
        private final double R = 6371000; // 距离计算使用的平均地球半径（米）；面积计算在 ringArea 使用 6378137

        // CRS缓存，避免重复解析WKT
        private final ConcurrentHashMap<String, CoordinateReferenceSystem> crsCache = new ConcurrentHashMap<>();
        // 变换缓存，避免重复构建同一分带的投影转换
        private final ConcurrentHashMap<String, MathTransform> txCache = new ConcurrentHashMap<>();

        // GeometryFactory缓存，避免重复创建
        private final GeometryFactory geometryFactory = new GeometryFactory();

        // 预定义的空几何体，避免重复创建
        private final Geometry EMPTY_GEOMETRY = geometryFactory.createGeometryCollection(null);

        // 空几何体的WKT表示，用于初始化空结果
        private final String EMPTY_WKT = "GEOMETRYCOLLECTION EMPTY";

        // 分批处理时的批次大小（点数量）
        private final int BATCH_SIZE = 500;

        // 默认轮廓返回的最多多边形数量（TopN）
        private final int DEFAULT_MAX_OUTLINE_SEGMENTS = 10;

        // 圆近似细分：让小面积轨迹缓冲区更圆滑（半圆效果）
        private final int DEFAULT_BUFFER_QUADRANT = 2;

        // 线缓冲样式（更圆滑边界）：拐角样式、端头样式与mitre限制
        private final int BUFFER_JOIN_STYLE = BufferParameters.JOIN_ROUND;
        // 线缓冲的端头样式（圆端），影响线段端点的缓冲形状
        private final int BUFFER_END_CAP_STYLE = BufferParameters.CAP_ROUND;
        // 线缓冲的mitre限制值，控制锐角的处理方式，值越小角越圆滑，值越大允许越尖的角
        private final double BUFFER_MITRE_LIMIT = 2.0;

        // 轨迹分段最小切割距离（米），用于splitRoad方法中的轨迹点分组
        private final TreeMap<Double, Double> MIN_SEGMENT_DISTANCE_THRESHOLD_MAP = new TreeMap<Double, Double>() {
            {
                put(1.0, 3.3);
                put(1.75, 2.7);
                put(2.8, 2.8);
                put(Double.MAX_VALUE, 4.0);
            }
        };

    }

    /**
     * 构建器类，用于构建GisUtil实例。
     * <p>
     * 实现构建器模式，允许灵活配置和创建GisUtil对象，提供流式API进行参数设置。
     * </p>
     */
    public static class Builder {
        // 配置对象，包含各种常量和默认值
        private final Config config = new Config();

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
     * 当前仅清理内部缓存；预留扩展以释放外部资源。
     */
    @Override
    public void close() {
        log.info("[销毁{}] 开始", this.getClass().getSimpleName());
        // 清理缓存
        config.crsCache.clear();
        // 回收各种资源（预留扩展点）
        log.info("[销毁{}] 结束", this.getClass().getSimpleName());
    }

    /**
     * 将实际时间间隔映射到标准频率（正负2秒容错）
     * 
     * @param actualInterval      实际间隔秒数
     * @param standardFrequencies 标准频率数组
     * @return 对应的标准频率，如果无法匹配返回-1
     */
    private static int mapToStandardFrequency(int actualInterval, int[] standardFrequencies) {
        for (int freq : standardFrequencies) {
            if (Math.abs(actualInterval - freq) <= 2) {
                return freq;
            }
        }
        return -1;
    }

    /**
     * 获取并缓存坐标参考系（CRS）。
     * 
     * <p>
     * 使用线程安全的并发缓存（ConcurrentHashMap）按坐标系代码惰性解析和存储CRS对象，
     * 避免重复调用昂贵的 {@code CRS.decode()} 操作，显著提升性能。
     * </p>
     * 
     * <p>
     * 支持的标准坐标系代码格式：
     * <ul>
     * <li>EPSG代码：EPSG:4326、EPSG:3857</li>
     * <li>ESRI代码：ESRI:54030</li>
     * <li>自定义WKT字符串</li>
     * </ul>
     * </p>
     * 
     * @param code 坐标系代码，如 "EPSG:4326"，不能为null或空字符串
     * @return 解析后的坐标参考系对象，保证线程安全
     * @throws IllegalArgumentException 如果code为null、空字符串或格式无效
     * @throws RuntimeException         如果CRS解析失败，包含详细的错误信息
     */
    private CoordinateReferenceSystem getCachedCRS(String code) {
        // 参数校验：确保坐标系代码有效
        if (code == null || code.trim().isEmpty()) {
            throw new IllegalArgumentException("坐标系代码不能为空");
        }

        // 使用computeIfAbsent确保线程安全的懒加载
        return config.crsCache.computeIfAbsent(code.trim(), key -> {
            try {
                // 使用lenient模式解析，提高兼容性
                return CRS.decode(key, true);
            } catch (Exception e) {
                // 提供更详细的错误信息，便于调试
                throw new RuntimeException(String.format("解析坐标参考系失败: code='%s', 错误: %s",
                        key, e.getMessage()), e);
            }
        });
    }

    /**
     * 按经度计算6度分带区号。
     * 规则：zone = floor(lon / 6) + 1，并限制下界为1（上界由使用场景保证）。
     *
     * @param lon 经度（度）
     * @return 分带区号
     */
    private int gkZoneFromLon(double lon) {
        // 6°分带：区号计算 floor(lon/6)+1；中央经线 L0 = 6*zone - 3
        int zone = (int) Math.floor(lon / 6.0) + 1;
        // 区号最小限制为 1（度制经度输入），不做国家范围特定限制
        if (zone < 1)
            zone = 1;
        return zone;
    }

    /**
     * 计算分带的中央经线（度）。
     * 公式：central_meridian = zone * 6 - 3。
     *
     * @param zone 分带区号
     * @return 中央经线（度）
     */
    private double gkCentralMeridian(int zone) {
        return zone * 6.0 - 3.0;
    }

    /**
     * 计算假东距（米）。
     * 约定：false_easting = zone * 1,000,000 + 500,000（仅影响坐标值，不影响几何形状）。
     *
     * @param zone 分带区号
     * @return 假东距（米）
     */
    private double gkFalseEasting(int zone) {
        // 通用约定：假东距 = 区号 * 1,000,000 + 500,000（不会影响形状，仅影响坐标值）
        return zone * 1_000_000.0 + 500_000.0;
    }

    /**
     * 按高斯投影东距近似反推分带区号。
     * 公式：round((E - 500000) / 1000000)；并限制在[1,60]。
     *
     * @param easting 东距（米）
     * @return 分带区号（1-60）
     */
    private int gkZoneFromEasting(double easting) {
        int zone = (int) Math.round((easting - 500_000.0) / 1_000_000.0);
        if (zone < 1)
            zone = 1;
        if (zone > 60)
            zone = 60;
        return zone;
    }

    /**
     * 从几何包络的质心东距推断分带。
     * 输入为空或包络无效时返回安全中间带30。
     *
     * @param g 输入几何（高斯坐标系假定）
     * @return 分带区号
     */
    private int inferGkZoneFromGeometryEasting(Geometry g) {
        if (g == null || g.isEmpty()) {
            return 30; // 安全回退：中间分带，避免异常
        }
        Envelope env = g.getEnvelopeInternal();
        if (env == null || env.isNull()) {
            return 30;
        }
        double easting = env.getMinX() + env.getWidth() / 2.0;
        return gkZoneFromEasting(easting);
    }

    /**
     * 构造并缓存高斯-克吕格（横轴墨卡托）投影CRS（6度分带）。
     * 基于WKT使用 WGS84 椭球，按分带设置中央经线与假东距，缓存键：GK6:<zone>。
     *
     * @param lon WGS84经度（度），用于选择分带
     * @return 对应分带的米制投影CRS
     */
    private CoordinateReferenceSystem getGaussKrugerCRSByLon(double lon) throws Exception {
        int zone = gkZoneFromLon(lon);
        double central = gkCentralMeridian(zone);
        double falseEasting = gkFalseEasting(zone);
        String key = "GK6:" + zone;
        CoordinateReferenceSystem crs = config.crsCache.get(key);
        if (crs != null)
            return crs;

        // 基于WKT定义高斯-克吕格（Transverse_Mercator）投影，基于WGS84椭球
        String wkt = "PROJCS[\"WGS 84 / Gauss-Kruger zone " + zone + "\", " +
                "GEOGCS[\"WGS 84\", DATUM[\"WGS_1984\", SPHEROID[\"WGS 84\",6378137,298.257223563]], " +
                "PRIMEM[\"Greenwich\",0], UNIT[\"degree\",0.0174532925199433]], " +
                "PROJECTION[\"Transverse_Mercator\"], " +
                "PARAMETER[\"latitude_of_origin\",0], " +
                "PARAMETER[\"central_meridian\"," + central + "], " +
                "PARAMETER[\"scale_factor\",1], " +
                "PARAMETER[\"false_easting\"," + falseEasting + "], " +
                "PARAMETER[\"false_northing\",0], " +
                "UNIT[\"metre\",1], AXIS[\"Easting\",EAST], AXIS[\"Northing\",NORTH]]";
        crs = CRS.parseWKT(wkt);
        config.crsCache.put(key, crs);
        return crs;
    }

    /**
     * 获取并缓存 WGS84→高斯（米制）前向变换。
     * 源CRS：EPSG:4326；目标CRS：按经度选择分带的高斯-克吕格；缓存键：W2G:<zone>。
     *
     * @param lon WGS84经度（度），用于选择分带
     * @return 前向投影变换（WGS84→Gauss）
     */
    private MathTransform getTxWgsToGaussByLon(double lon) throws Exception {
        int zone = gkZoneFromLon(lon);
        String key = "W2G:" + zone;
        MathTransform cached = config.txCache.get(key);
        if (cached != null)
            return cached;
        CoordinateReferenceSystem src = getCachedCRS(config.WGS84);
        CoordinateReferenceSystem tgt = getGaussKrugerCRSByLon(lon);
        MathTransform mt = CRS.findMathTransform(src, tgt, true);
        config.txCache.put(key, mt);
        return mt;
    }

    /**
     * 获取并缓存 高斯→WGS84 逆向变换。
     * 源CRS：按经度选择分带的高斯-克吕格；目标CRS：EPSG:4326；缓存键：G2W:<zone>。
     *
     * @param wgs84Lon WGS84经度（度），用于选择分带
     * @return 逆向投影变换（Gauss→WGS84）
     */
    private MathTransform getTxGaussToWgsByLon(double wgs84Lon) throws Exception {
        int zone = gkZoneFromLon(wgs84Lon);
        String key = "G2W:" + zone;
        MathTransform cached = config.txCache.get(key);
        if (cached != null)
            return cached;
        CoordinateReferenceSystem src = getGaussKrugerCRSByLon(wgs84Lon);
        CoordinateReferenceSystem tgt = getCachedCRS(config.WGS84);
        MathTransform mt = CRS.findMathTransform(src, tgt, true);
        config.txCache.put(key, mt);
        return mt;
    }

    /**
     * 校验 WGS84 坐标几何的经纬度范围是否合理。
     * 判断规则：经度在 [-180, 180]、纬度在 [-90, 90]，且坐标不为 (0,0)；不显式检查 NaN/Inf（解析阶段应已过滤）。
     *
     * @param g WGS84 坐标系下的几何
     * @return 坐标范围合理返回 true，否则 false
     */
    private boolean isValidWgs84Geometry(Geometry g) {
        Coordinate[] coordinates = g.getCoordinates();
        for (Coordinate coord : coordinates) {
            // 检查经度是否在[-180, 180]范围内
            if (Math.abs(coord.x) > 180) {
                return false;
            }
            // 检查纬度是否在[-90, 90]范围内
            if (Math.abs(coord.y) > 90) {
                return false;
            }
            // 检查坐标是否为0,0（通常为无效坐标）
            if (coord.x == 0 && coord.y == 0) {
                return false;
            }
        }
        return true;
    }

    /**
     * 过滤异常点并按时间升序排序。
     * 条件：时间非空；经纬度在有效范围；排除(0,0)；保留 `isValidWgs84Geometry` 检查通过的点。
     *
     * @param seg 原始轨迹点列表
     * @return 过滤并按时间排序后的轨迹点列表
     */
    private List<TrackPoint> filterAndSortTrackPoints(List<TrackPoint> seg) {
        if (seg == null)
            return Collections.emptyList();
        Comparator<LocalDateTime> cmp = Comparator
                .nullsLast(Comparator.naturalOrder());
        return seg.stream()
                .filter(p -> p.getTime() != null)
                .filter(p -> Math.abs(p.getLon()) <= 180 && Math.abs(p.getLat()) <= 90)
                .filter(p -> !(p.getLon() == 0 && p.getLat() == 0))
                .sorted(Comparator.comparing(TrackPoint::getTime, cmp))
                .collect(Collectors.toList());
    }

    /**
     * 检查几何坐标的基本有效性。
     * 判断规则：坐标不为 (0,0)，且不存在接近 0 的可疑值（|x|<1e-10 或 |y|<1e-10）；不负责空几何与 NaN/Inf
     * 检查（调用方需在前置流程处理）。
     *
     * @param g 几何（应为非空）
     * @return 坐标通过基本检查返回 true，否则 false
     */
    private boolean hasValidCoordinates(Geometry g) {
        Coordinate[] coordinates = g.getCoordinates();
        for (Coordinate coord : coordinates) {
            // 检查是否存在非常接近0但非0的坐标值，这通常表示转换错误
            if ((Math.abs(coord.x) < 1e-10 && Math.abs(coord.x) > 0) ||
                    (Math.abs(coord.y) < 1e-10 && Math.abs(coord.y) > 0)) {
                return false;
            }
            // 检查坐标是否为0,0（通常为无效坐标）
            if (coord.x == 0 && coord.y == 0) {
                return false;
            }
        }
        return true;
    }

    /**
     * 计算球面环面积（平方米），与 Turf.js 的 ringArea 对齐。
     * 使用 WGS84 半径（6378137）与经纬度弧度。
     *
     * @param ring 外/内环线
     * @return 面积（平方米）
     */
    private static double ringArea(LineString ring) {
        Coordinate[] coords = ring.getCoordinates();
        int len = (coords == null) ? 0 : coords.length;
        if (len <= 2)
            return 0.0;
        double area = 0.0;
        for (int i = 0; i < len; i++) {
            Coordinate p1, p2, p3;
            if (i == len - 2) {
                p1 = coords[i];
                p2 = coords[i + 1];
                p3 = coords[0];
            } else if (i == len - 1) {
                p1 = coords[i];
                p2 = coords[0];
                p3 = coords[1];
            } else {
                p1 = coords[i];
                p2 = coords[i + 1];
                p3 = coords[i + 2];
            }
            area += (toRad(p3.x) - toRad(p1.x)) * Math.sin(toRad(p2.y));
        }
        double R = 6378137.0; // Turf 默认采用的半径（WGS84）
        return area * R * R / 2.0;
    }

    /**
     * 角度转弧度。
     *
     * @param deg 角度
     * @return 弧度
     */
    private static double toRad(double deg) {
        return deg * Math.PI / 180.0;
    }

    /**
     * 简化版线缓冲区构建方法（无分段逻辑）。
     * 
     * <p>
     * 与 {@link #buildSmartLineBuffer} 方法相比，此方法功能简化：
     * - 不做轨迹分段判断，直接处理全部轨迹点
     * - 不做道格拉斯-普克简化，保留原始轨迹精度
     * - 适用于轨迹较短或无需分段处理的场景
     * </p>
     * 
     * <p>
     * 处理流程：
     * 1. WGS84坐标转高斯-克吕格投影
     * 2. 构建线几何（保留所有点）
     * 3. 创建指定宽度的缓冲区
     * 4. 转换回WGS84坐标系
     * </p>
     * 
     * @param seg    轨迹点列表，每个点包含经纬度坐标
     * @param widthM 缓冲区宽度，单位为米
     * @return 缓冲区几何图形（WGS84坐标系）
     * @throws Exception 当坐标转换或几何处理失败时抛出
     * 
     * @see #buildSmartLineBuffer(List, double, boolean) 智能版缓冲区构建方法（支持分段）
     */
    private Geometry buildSimpleLineBuffer(List<TrackPoint> seg, double widthM) throws Exception {
        // 参数校验
        if (seg == null || seg.size() < 2) {
            return config.geometryFactory.createGeometryCollection(null);
        }

        long startTime = System.currentTimeMillis();
        log.trace("[buildSmartLineBuffer] 开始处理 {} 个轨迹点，线缓冲宽度: {} 米", seg.size(), widthM);

        // 获取原点经度并创建坐标转换
        double originLon = seg.get(0).getLon();
        MathTransform txWgsToGk = getTxWgsToGaussByLon(originLon);

        // 将轨迹点转换为高斯投影坐标
        List<Coordinate> coordinates = new ArrayList<>(seg.size());
        for (TrackPoint point : seg) {
            Coordinate src = new Coordinate(point.getLon(), point.getLat());
            Coordinate target = new Coordinate();
            JTS.transform(src, target, txWgsToGk);
            coordinates.add(target);
        }

        // 创建线几何（保留所有轨迹点，不做简化）
        Coordinate[] coordArray = coordinates.toArray(new Coordinate[0]);
        LineString line = config.geometryFactory.createLineString(coordArray);

        // 设置缓冲区参数
        BufferParameters params = new BufferParameters();
        params.setQuadrantSegments(config.DEFAULT_BUFFER_QUADRANT);
        params.setEndCapStyle(config.BUFFER_END_CAP_STYLE);
        params.setJoinStyle(config.BUFFER_JOIN_STYLE);
        params.setMitreLimit(config.BUFFER_MITRE_LIMIT);

        // 创建缓冲区（使用原始轨迹点，保证精度）
        Geometry buffered = BufferOp.bufferOp(line, widthM, params);

        // 转换回WGS84坐标系
        MathTransform txGkToWgs = getTxGaussToWgsByLon(originLon);
        Geometry result = JTS.transform(buffered, txGkToWgs);

        // 记录处理时间
        long endTime = System.currentTimeMillis();
        log.debug("[buildSimpleLineBuffer] 完成，总计耗时: {}ms", endTime - startTime);

        return result;
    }

    /**
     * WGS84坐标系转换为高斯投影坐标系
     * 
     * @param wgs84Points WGS84坐标系的轨迹点列表
     * @return 高斯投影坐标系的轨迹点列表
     */
    private List<TrackPoint> convertWgs84ToGaussProjection(List<TrackPoint> wgs84Points) {
        if (CollUtil.isEmpty(wgs84Points)) {
            return new ArrayList<>();
        }

        List<TrackPoint> result = new ArrayList<>();

        try {
            // 获取高斯投影坐标系
            CoordinateReferenceSystem gaussCRS = getGaussKrugerCRSByLon(wgs84Points.get(0).getLon());
            CoordinateReferenceSystem wgs84CRS = DefaultGeographicCRS.WGS84;

            // 创建坐标转换
            MathTransform transform = CRS.findMathTransform(wgs84CRS, gaussCRS, true);

            for (TrackPoint point : wgs84Points) {
                try {
                    // 创建WGS84坐标
                    DirectPosition2D wgs84Pos = new DirectPosition2D(wgs84CRS, point.getLon(), point.getLat());

                    // 转换为高斯投影坐标
                    DirectPosition2D gaussPos = new DirectPosition2D();
                    transform.transform(wgs84Pos, gaussPos);

                    // 创建新的轨迹点，使用高斯投影坐标（单位：米）
                    TrackPoint gaussPoint = new TrackPoint(point.getTime(), gaussPos.getX(), gaussPos.getY());

                    result.add(gaussPoint);

                } catch (Exception e) {
                    log.warn("单点坐标转换失败: lat={}, lon={}, error: {}",
                            point.getLat(), point.getLon(), e.getMessage());
                    // 跳过转换失败的点
                }
            }

            log.debug("WGS84转高斯投影完成: 原始点数={}, 转换成功点数={}",
                    wgs84Points.size(), result.size());

        } catch (Exception e) {
            log.warn("坐标转换初始化失败", e.getMessage());
            return new ArrayList<>();
        }

        return result;
    }

    /**
     * 在高斯投影坐标系下，过滤出在给定几何图形内的轨迹点
     * 渐进式精确过滤：边界框预过滤 → 多线程并行处理 → 100%精确验证
     * 支持多轮廓、复杂几何、确保零误判
     * 
     * @param gaussPoints 高斯投影坐标系的轨迹点列表（单位：米）
     * @param filterGeom  高斯投影坐标系下的过滤几何图形（支持多轮廓）
     * @return 高斯投影坐标系的轨迹点列表（仅包含在filterGeom内的点）
     */
    private List<TrackPoint> filterGaussPointsByGeometry(List<TrackPoint> gaussPoints, Geometry filterGeom) {
        if (CollUtil.isEmpty(gaussPoints) || filterGeom == null || filterGeom.isEmpty()) {
            return new ArrayList<>();
        }

        // ===== 第一阶段：边界框预过滤（内存高效）=====
        Envelope filterEnv = filterGeom.getEnvelopeInternal();
        List<TrackPoint> bboxFiltered = new ArrayList<>();

        for (TrackPoint point : gaussPoints) {
            if (filterEnv.contains(point.getLon(), point.getLat())) {
                bboxFiltered.add(point);
            }
        }

        log.debug("边界框预过滤: 原始点数={}, 边界框内点数={}, 过滤几何类型={}",
                gaussPoints.size(), bboxFiltered.size(), filterGeom.getGeometryType());

        // 如果边界框内无点，直接返回
        if (bboxFiltered.isEmpty()) {
            return bboxFiltered;
        }

        // ===== 第二阶段：多线程精确过滤（确保100%准确）=====
        return preciseParallelFilter(bboxFiltered, filterGeom);
    }

    /**
     * 多线程精确空间过滤（100%准确，零误判）
     * 支持复杂几何、多轮廓、确保每个点都经过精确几何判断
     */
    private List<TrackPoint> preciseParallelFilter(List<TrackPoint> points, Geometry filterGeom) {
        // 小数据集：单线程处理（避免线程开销）
        if (points.size() < 1000) {
            return preciseSingleThreadFilter(points, filterGeom);
        }

        // 大数据集：多线程批处理
        return preciseMultiThreadFilter(points, filterGeom);
    }

    /**
     * 单线程精确过滤（小数据集）
     */
    private List<TrackPoint> preciseSingleThreadFilter(List<TrackPoint> points, Geometry filterGeom) {
        List<TrackPoint> result = new ArrayList<>();

        for (TrackPoint point : points) {
            try {
                Point pointGeom = config.geometryFactory.createPoint(
                        new Coordinate(point.getLon(), point.getLat()));

                // 精确几何判断：支持复杂多轮廓、孔洞等
                if (filterGeom.contains(pointGeom)) {
                    result.add(point);
                }
            } catch (Exception e) {
                log.warn("单线程空间过滤失败: x={}, y={}, error: {}",
                        point.getLon(), point.getLat(), e.getMessage());
            }
        }

        log.trace("单线程精确过滤完成: 输入点数={}, 结果点数={}", points.size(), result.size());
        return result;
    }

    /**
     * 多线程精确过滤（大数据集）
     * 批处理 + 并行流，确保高性能和准确性
     */
    private List<TrackPoint> preciseMultiThreadFilter(List<TrackPoint> points, Geometry filterGeom) {
        final int THREAD_BATCH_SIZE = Math.max(100, points.size() / (Runtime.getRuntime().availableProcessors() * 4));

        List<TrackPoint> result = Collections.synchronizedList(new ArrayList<>());

        // 按批次并行处理
        IntStream.range(0, (points.size() + THREAD_BATCH_SIZE - 1) / THREAD_BATCH_SIZE)
                .parallel()
                .forEach(batchIndex -> {
                    int start = batchIndex * THREAD_BATCH_SIZE;
                    int end = Math.min(start + THREAD_BATCH_SIZE, points.size());
                    List<TrackPoint> batch = points.subList(start, end);

                    // 处理单个批次
                    List<TrackPoint> batchResult = preciseSingleThreadFilter(batch, filterGeom);
                    result.addAll(batchResult);

                    log.trace("批次处理完成: 批次索引={}, 批次大小={}, 结果大小={}",
                            batchIndex, batch.size(), batchResult.size());
                });

        log.debug("多线程精确过滤完成: 输入点数={}, 结果点数={}, 批次大小={}, 处理器数={}",
                points.size(), result.size(), THREAD_BATCH_SIZE, Runtime.getRuntime().availableProcessors());

        return new ArrayList<>(result);
    }

    /**
     * 使用WGS84坐标系的WKT字符串计算面积（亩）
     * 采用与Turf.js相同的球面面积计算算法，确保计算结果一致
     * 
     * @param wkt WGS84坐标系的WKT字符串
     * @return 面积（亩），四舍五入保留4位小数；非法或空输入返回0
     */
    private double calcMuByWgs84Wkt(String wkt) {
        if (StrUtil.isBlank(wkt)) {
            return 0.0;
        }

        try {
            // 解析WKT字符串为几何图形
            WKTReader reader = new WKTReader(config.geometryFactory);
            Geometry geometry = reader.read(wkt);

            if (geometry == null || geometry.isEmpty()) {
                return 0.0;
            }

            // 使用GeoTools的球面面积计算，与Turf.js算法一致
            double areaSqm = calculateSphericalArea(geometry);

            // 转换为亩：1亩 = 2000/3平方米，四舍五入保留4位小数
            return Math.round((areaSqm / (2000.0 / 3.0)) * 10000.0) / 10000.0;

        } catch (Exception e) {
            log.warn("WKT字符串计算面积失败: {}", e.getMessage());
            return 0.0;
        }
    }

    /**
     * 计算球面面积（平方米）
     * 使用GeoTools的球面面积计算算法，与Turf.js的球面面积计算结果一致
     * 
     * @param geometry WGS84坐标系的几何图形
     * @return 球面面积（平方米）
     */
    private double calculateSphericalArea(Geometry geometry) {
        if (geometry == null || geometry.isEmpty()) {
            return 0.0;
        }

        double totalArea = 0.0;

        if (geometry instanceof Polygon) {
            totalArea = calculatePolygonSphericalArea((Polygon) geometry);
        } else if (geometry instanceof MultiPolygon) {
            MultiPolygon multiPolygon = (MultiPolygon) geometry;
            for (int i = 0; i < multiPolygon.getNumGeometries(); i++) {
                Polygon polygon = (Polygon) multiPolygon.getGeometryN(i);
                totalArea += calculatePolygonSphericalArea(polygon);
            }
        }

        return Math.abs(totalArea); // 确保面积为正值
    }

    /**
     * 计算单个多边形的球面面积（平方米）
     * 使用球面多边形面积公式，考虑地球曲率
     * 
     * @param polygon WGS84坐标系的多边形
     * @return 球面面积（平方米）
     */
    private double calculatePolygonSphericalArea(Polygon polygon) {
        // 外环面积
        double exteriorArea = calculateRingSphericalArea(polygon.getExteriorRing());

        // 减去内环（孔洞）面积
        double holesArea = 0.0;
        for (int i = 0; i < polygon.getNumInteriorRing(); i++) {
            holesArea += calculateRingSphericalArea(polygon.getInteriorRingN(i));
        }

        return exteriorArea - holesArea;
    }

    /**
     * 计算线环的球面面积（平方米）
     * 使用球面多边形面积公式：A = R² * |Σ(λi+1 - λi) * sin(φi+1 + φi)/2|
     * 其中R为地球半径，λ为经度，φ为纬度
     * 
     * @param ring 线环（WGS84坐标系，单位：度）
     * @return 球面面积（平方米）
     */
    private double calculateRingSphericalArea(LineString ring) {
        if (ring == null || ring.isEmpty()) {
            return 0.0;
        }

        // WGS84椭球长半轴（米），与Turf.js保持一致
        final double EARTH_RADIUS = 6378137.0;

        Coordinate[] coords = ring.getCoordinates();
        if (coords.length < 3) {
            return 0.0; // 需要至少3个点才能形成多边形
        }

        double area = 0.0;

        // 使用球面多边形面积公式
        for (int i = 0; i < coords.length - 1; i++) {
            double lon1 = Math.toRadians(coords[i].x);
            double lat1 = Math.toRadians(coords[i].y);
            double lon2 = Math.toRadians(coords[i + 1].x);
            double lat2 = Math.toRadians(coords[i + 1].y);

            // 球面多边形面积公式
            area += (lon2 - lon1) * (2 + Math.sin(lat1) + Math.sin(lat2));
        }

        area = Math.abs(area) * EARTH_RADIUS * EARTH_RADIUS / 2.0;
        return area;
    }

    /**
     * 将高斯投影坐标系的轨迹点列表转换为WGS84坐标系
     * 
     * @param gaussPoints 高斯投影坐标系的轨迹点列表（单位：米）
     * @param wgs84Lon    WGS84经度（度），用于确定高斯投影带；通常使用轨迹点的经度
     * @return WGS84坐标系的轨迹点列表；转换失败时返回空列表
     */
    private List<TrackPoint> gaussPointsToWgs84(List<TrackPoint> gaussPoints, double wgs84Lon) {
        if (CollUtil.isEmpty(gaussPoints)) {
            return new ArrayList<>();
        }

        List<TrackPoint> result = new ArrayList<>();

        try {
            // 获取高斯投影坐标系
            CoordinateReferenceSystem gaussCRS = getGaussKrugerCRSByLon(wgs84Lon);
            CoordinateReferenceSystem wgs84CRS = DefaultGeographicCRS.WGS84;

            // 创建坐标转换（高斯投影到WGS84）
            MathTransform transform = CRS.findMathTransform(gaussCRS, wgs84CRS, true);

            for (TrackPoint gaussPoint : gaussPoints) {
                try {
                    // 创建高斯投影坐标
                    DirectPosition2D gaussPos = new DirectPosition2D(gaussCRS, gaussPoint.getLon(),
                            gaussPoint.getLat());

                    // 转换为WGS84坐标
                    DirectPosition2D wgs84Pos = new DirectPosition2D();
                    transform.transform(gaussPos, wgs84Pos);

                    // 创建新的轨迹点，使用WGS84坐标（单位：度）
                    TrackPoint wgs84Point = new TrackPoint(gaussPoint.getTime(), wgs84Pos.getX(), wgs84Pos.getY());

                    result.add(wgs84Point);

                } catch (Exception e) {
                    log.warn("单点坐标转换失败（高斯投影→WGS84）: x={}, y={}, error: {}",
                            gaussPoint.getLon(), gaussPoint.getLat(), e.getMessage());
                    // 跳过转换失败的点
                }
            }

            log.debug("高斯投影转WGS84完成: 原始点数={}, 转换成功点数={}",
                    gaussPoints.size(), result.size());

        } catch (Exception e) {
            log.warn("高斯投影转WGS84坐标转换初始化失败", e);
            return new ArrayList<>();
        }

        return result;
    }

    /**
     * 将高斯投影坐标系的几何图形转换为WGS84坐标系的WKT字符串
     * 
     * @param gaussGeometry 高斯投影坐标系的几何图形
     * @param wgs84Lon      WGS84经度（度），用于确定高斯投影带；通常使用轨迹点的经度
     * @return WGS84坐标系的WKT字符串；转换失败时返回空几何的WKT
     */
    private String gaussGeometryToWgs84Wkt(Geometry gaussGeometry, double wgs84Lon) {
        if (gaussGeometry == null || gaussGeometry.isEmpty()) {
            return config.EMPTY_WKT;
        }

        try {
            // 获取高斯投影到WGS84的坐标转换
            MathTransform txBack = getTxGaussToWgsByLon(wgs84Lon);

            // 将高斯投影几何图形转换到WGS84坐标系
            Geometry wgs84Geometry = JTS.transform(gaussGeometry, txBack);

            // 转换为WKT字符串
            return wgs84Geometry.toText();
        } catch (Exception e) {
            log.warn("高斯投影几何图形转WGS84 WKT失败: {}", e.getMessage());
            return config.EMPTY_WKT;
        }
    }

    /**
     * 计算两点球面距离（哈弗辛公式，WGS84 平均地球半径）。
     * 输入为经纬度（度），返回沿地球表面的近似测地距离（米）。
     * 注意：这是球面近似的地理距离，不是投影平面直线距离；若需平面距离，请先做米制投影后再计算。
     *
     * @param p1 点1（WGS84，经纬度，单位度）
     * @param p2 点2（WGS84，经纬度，单位度）
     * @return 距离（米）
     */
    public double haversine(CoordinatePoint p1, CoordinatePoint p2) {
        double dLat = Math.toRadians(p2.getLat() - p1.getLat());
        double dLon = Math.toRadians(p2.getLon() - p1.getLon());
        double a = Math.sin(dLat / 2) * Math.sin(dLat / 2) +
                Math.cos(Math.toRadians(p1.getLat())) *
                        Math.cos(Math.toRadians(p2.getLat())) *
                        Math.sin(dLon / 2) * Math.sin(dLon / 2);
        return config.R * 2 * Math.atan2(Math.sqrt(a), Math.sqrt(1 - a));
    }

    /**
     * 点是否在多边形内（WGS84）
     * 接收点（WGS84，经纬度）和 WKT（WGS84，`POLYGON`/`MULTIPOLYGON`），判断点是否在多边形内；边界视为内。
     *
     * @param point      点坐标（WGS84，经度 `lon`、纬度 `lat`，单位度）
     * @param wktPolygon 多边形 WKT（WGS84，类型为 `POLYGON` 或 `MULTIPOLYGON`）
     * @return 是否在内（含边界）；解析失败、类型不支持或为空几何返回 false
     *
     * @implNote 先用包络矩形快速裁剪，再用 `PreparedGeometry#covers(Point)` 判断；若多边形坐标非法，先通过
     *           `buffer(0)` 进行拓扑修复。
     */
    public boolean pointInPolygon(CoordinatePoint point, String wktPolygon) {
        if (point == null || wktPolygon == null || wktPolygon.trim().isEmpty()) {
            return false;
        }
        double lon = point.getLon();
        double lat = point.getLat();
        if (Double.isNaN(lon) || Double.isNaN(lat)) {
            return false;
        }
        if (lon < -180.0 || lon > 180.0 || lat < -90.0 || lat > 90.0) {
            return false;
        }
        try {
            WKTReader reader = new WKTReader(config.geometryFactory);
            Geometry poly = reader.read(wktPolygon);
            if (poly == null || poly.isEmpty()) {
                return false;
            }
            String type = poly.getGeometryType().toUpperCase();
            if (!("POLYGON".equals(type) || "MULTIPOLYGON".equals(type))) {
                return false;
            }
            if (!hasValidCoordinates(poly)) {
                poly = poly.buffer(0);
                if (poly.isEmpty()) {
                    return false;
                }
            }
            Geometry pt = config.geometryFactory.createPoint(new Coordinate(lon, lat));
            Envelope env = poly.getEnvelopeInternal();
            if (!env.contains(pt.getCoordinate())) {
                return false;
            }
            PreparedGeometry prep = PreparedGeometryFactory.prepare(poly);
            return prep.covers(pt);
        } catch (ParseException ex) {
            log.error("解析 WKT 失败: {}", ex.getMessage());
            return false;
        } catch (Exception ex) {
            log.error("pointInPolygon 判断出错: {}", ex.getMessage(), ex);
            return false;
        }
    }

    /**
     * 将几何转换为WKT（统一WGS84）
     * 自动识别坐标系：若坐标范围合理则视为 WGS84 直接输出；否则按高斯-克吕格投影推断分带并转换到 WGS84 后输出。
     * 处理策略：null 返回 GEOMETRYCOLLECTION EMPTY；空几何直接 `toText()`；坐标非法时先 `buffer(0)`
     * 修复再输出。
     *
     * @param g 几何对象（WGS84 或高斯-克吕格），支持任意 `Geometry`
     * @return 统一为 WGS84 的 WKT 字符串
     *
     * @implNote 分带依据几何的东距范围推断，中央子午线由分带计算；异常时回退为原始 WKT 或空集合。
     */
    public String toWkt(Geometry g) {
        try {
            // 空几何直接返回其WKT
            if (g == null) {
                return config.EMPTY_WKT;
            }
            if (g.isEmpty()) {
                return g.toText();
            }
            // 自动识别坐标系：WGS84直接输出，否则视为高斯并转换到WGS84
            Geometry wgs;
            if (isValidWgs84Geometry(g)) {
                wgs = g;
            } else {
                int zone = inferGkZoneFromGeometryEasting(g);
                double lon = gkCentralMeridian(zone);
                MathTransform tx = getTxGaussToWgsByLon(lon);
                wgs = JTS.transform(g, tx);
            }
            if (!hasValidCoordinates(wgs)) {
                wgs = wgs.buffer(0);
            }
            return wgs.toText();
        } catch (Exception e) {
            log.error("坐标转换到WGS84失败: " + e.getMessage(), e);
            return g != null ? g.toText() : config.EMPTY_WKT;
        }
    }

    /**
     * 计算几何图形的面积（mu单位，支持 Polygon 与 MultiPolygon）。
     * - 自动识别坐标系：若为 WGS84（经纬度范围合理）则直接按球面公式计算；否则视为高斯-克吕格投影并转换到 WGS84。
     * - 面积计算与 Turf.js 对齐（球面面积，平方米）。
     *
     * @param outline 几何图形（POLYGON 或 MULTIPOLYGON）
     *
     * @return 面积（mu），以亩为单位，保留4位小数
     *
     * @throws RuntimeException 如果面积计算过程中发生错误
     */
    public double calcMu(Geometry outline) throws RuntimeException {
        if (outline == null || outline.isEmpty()) {
            return 0.0;
        }
        try {
            // 自动识别坐标系：WGS84 直接使用；否则按高斯→WGS84转换
            Geometry wgs;
            if (isValidWgs84Geometry(outline)) {
                wgs = outline;
            } else {
                int zone = inferGkZoneFromGeometryEasting(outline);
                double lon = gkCentralMeridian(zone);
                MathTransform tx = getTxGaussToWgsByLon(lon);
                wgs = JTS.transform(outline, tx);
            }

            double areaSqm = 0.0;
            if (wgs instanceof Polygon) {
                Polygon p = (Polygon) wgs;
                double areaOuter = Math.abs(ringArea(p.getExteriorRing()));
                double holesArea = 0.0;
                for (int i = 0; i < p.getNumInteriorRing(); i++) {
                    holesArea += Math.abs(ringArea(p.getInteriorRingN(i)));
                }
                areaSqm = areaOuter - holesArea;
            } else if (wgs instanceof MultiPolygon) {
                MultiPolygon mp = (MultiPolygon) wgs;
                for (int i = 0; i < mp.getNumGeometries(); i++) {
                    Polygon p = (Polygon) mp.getGeometryN(i);
                    double areaOuter = Math.abs(ringArea(p.getExteriorRing()));
                    double holesArea = 0.0;
                    for (int j = 0; j < p.getNumInteriorRing(); j++) {
                        holesArea += Math.abs(ringArea(p.getInteriorRingN(j)));
                    }
                    areaSqm += (areaOuter - holesArea);
                }
            } else {
                return 0.0;
            }
            // 1 亩 = 666.6667 平方米
            return Math.round((areaSqm / 666.6667) * 10000.0) / 10000.0;
        } catch (Exception e) {
            throw new RuntimeException("计算面积时出错: " + e.getMessage(), e);
        }
    }

    /**
     * WKT面积计算（亩，WGS84）
     * 解析 `POLYGON`/`MULTIPOLYGON` 的 WKT（需为 WGS84），按球面公式计算面积并换算为亩，与 Turf.js 口径对齐。
     *
     * @param wkt WKT 字符串（WGS84），类型必须是 `POLYGON` 或 `MULTIPOLYGON`
     * @return 面积（亩），四舍五入保留 4 位小数；非法或空输入返回 0
     * @throws RuntimeException 解析失败或 WKT 非法时抛出
     *
     * @implNote 面积以外环减内环之和计算；环面积采用 WGS84 半径与经纬度弧度（参见 `ringArea`）。
     */
    public double calcMu(String wkt) {
        if (wkt == null || wkt.trim().isEmpty()) {
            return 0.0;
        }
        try {
            WKTReader reader = new WKTReader(config.geometryFactory);
            Geometry wgs = reader.read(wkt);
            double areaSqm = 0.0;
            if (wgs instanceof Polygon) {
                Polygon p = (Polygon) wgs;
                double areaOuter = Math.abs(ringArea(p.getExteriorRing()));
                double holesArea = 0.0;
                for (int i = 0; i < p.getNumInteriorRing(); i++) {
                    holesArea += Math.abs(ringArea(p.getInteriorRingN(i)));
                }
                areaSqm = areaOuter - holesArea;
            } else if (wgs instanceof MultiPolygon) {
                MultiPolygon mp = (MultiPolygon) wgs;
                for (int i = 0; i < mp.getNumGeometries(); i++) {
                    Polygon p = (Polygon) mp.getGeometryN(i);
                    double areaOuter = Math.abs(ringArea(p.getExteriorRing()));
                    double holesArea = 0.0;
                    for (int j = 0; j < p.getNumInteriorRing(); j++) {
                        holesArea += Math.abs(ringArea(p.getInteriorRingN(j)));
                    }
                    areaSqm += (areaOuter - holesArea);
                }
            } else {
                return 0.0;
            }
            return Math.round((areaSqm / 666.6667) * 10000.0) / 10000.0;
        } catch (ParseException e) {
            throw new RuntimeException("解析WKT字符串时出错: " + e.getMessage(), e);
        }
    }

    /**
     * 解析 WKT（WGS84）为 Geometry，并转换到高斯-克吕格米制坐标
     * 支持 `POLYGON` 与 `MULTIPOLYGON`；用于后续形态学运算与裁剪。
     *
     * @param wkt WKT 字符串（WGS84），类型必须为 `POLYGON` 或 `MULTIPOLYGON`
     * @return 高斯-克吕格投影下的几何对象（中央子午线依据几何质心经度选择）
     * @throws IllegalArgumentException WKT 为空或类型不支持时抛出
     * @throws RuntimeException         解析失败或坐标转换异常时抛出
     *
     * @implNote 投影转换缓存复用，避免频繁构建；转换失败时抛异常。
     */
    public Geometry fromWkt(String wkt) {
        if (wkt == null || wkt.trim().isEmpty()) {
            throw new IllegalArgumentException("WKT字符串不能为空");
        }
        try {
            WKTReader wktReader = new WKTReader(config.geometryFactory);
            Geometry wgs = wktReader.read(wkt);
            if (wgs == null) {
                throw new RuntimeException("无法解析WKT字符串: " + wkt);
            }
            String geometryType = wgs.getGeometryType().toUpperCase();
            if (!("POLYGON".equals(geometryType) || "MULTIPOLYGON".equals(geometryType))) {
                throw new IllegalArgumentException("不支持的几何类型: " + geometryType + "。仅支持POLYGON和MULTIPOLYGON");
            }
            double lon = wgs.getCentroid().getX();
            MathTransform tx = getTxWgsToGaussByLon(lon);
            Geometry gauss = JTS.transform(wgs, tx);
            return gauss;
        } catch (ParseException e) {
            throw new RuntimeException("解析WKT字符串时出错: " + e.getMessage(), e);
        } catch (Exception e) {
            throw new RuntimeException("WGS84→Gauss投影转换失败: " + e.getMessage(), e);
        }
    }

    /**
     * WKT 相交（WGS84）
     * 计算两个 WKT 的相交部分，返回相交 WKT（统一为 WGS84）与面积（亩）。
     *
     * @param wktA 输入几何 A 的 WKT（WGS84，`POLYGON`/`MULTIPOLYGON`）
     * @param wktB 输入几何 B 的 WKT（WGS84，`POLYGON`/`MULTIPOLYGON`）
     * @return `WktIntersectionResult`：当无相交或为空几何时 `wkt=null, mu=0`
     * @throws Exception 解析或坐标转换、几何运算失败时抛出
     *
     * @implNote 内部将 WKT 转为高斯-克吕格进行运算，再统一输出 WGS84；面积计算对齐球面公式。
     */
    public WktIntersectionResult intersection(String wktA, String wktB) throws Exception {
        Geometry g1 = fromWkt(wktA);
        Geometry g2 = fromWkt(wktB);

        // 坐标/几何修复（自相交等）
        if (!hasValidCoordinates(g1)) {
            g1 = g1.buffer(0);
            if (g1.isEmpty()) {
                return new WktIntersectionResult(null, 0.0);
            }
        }
        if (!hasValidCoordinates(g2)) {
            g2 = g2.buffer(0);
            if (g2.isEmpty()) {
                return new WktIntersectionResult(null, 0.0);
            }
        }

        // 几何简化以减少精度问题
        double tolerance = 1e-8; // WGS84坐标系下的简化容差
        g1 = DouglasPeuckerSimplifier.simplify(g1, tolerance);
        g2 = DouglasPeuckerSimplifier.simplify(g2, tolerance);

        Geometry inter;
        try {
            // 使用SnapIfNeededOverlayOp处理拓扑异常
            SnapIfNeededOverlayOp snapOp = new SnapIfNeededOverlayOp(g1, g2);
            inter = snapOp.getResultGeometry(OverlayOp.INTERSECTION);
        } catch (TopologyException e) {
            log.warn("SnapIfNeededOverlayOp失败，尝试使用缓冲交集: {}", e.getMessage());
            // 如果SnapIfNeededOverlayOp失败，使用缓冲交集作为备选方案
            double bufferDistance = 1e-10; // 极小的缓冲距离
            Geometry g1Buffer = g1.buffer(bufferDistance);
            Geometry g2Buffer = g2.buffer(bufferDistance);
            inter = g1Buffer.intersection(g2Buffer);
        }

        if (inter == null || inter.isEmpty()) {
            return new WktIntersectionResult(null, 0.0);
        }
        String wkt = toWkt(inter);
        double mu = calcMuByWgs84Wkt(wkt);
        return new WktIntersectionResult(wkt, mu);
    }

    /**
     * 判断两个 WKT（WGS84）是否相交
     * 接收两个 WKT 字符串（WGS84），解析为 `POLYGON`/`MULTIPOLYGON` 后判断是否相交。
     *
     * @param wktA WKT（WGS84），类型必须为 `POLYGON` 或 `MULTIPOLYGON`
     * @param wktB WKT（WGS84），类型必须为 `POLYGON` 或 `MULTIPOLYGON`
     * @return 是否相交；解析失败、类型不支持或为空几何返回 false
     *
     * @implNote 优先用包络矩形快速裁剪；必要时使用 PreparedGeometry 提升判断性能。
     */
    public boolean intersects(String wktA, String wktB) {
        if (wktA == null || wktA.trim().isEmpty() || wktB == null || wktB.trim().isEmpty()) {
            return false;
        }
        try {
            WKTReader reader = new WKTReader(config.geometryFactory);
            Geometry g1 = reader.read(wktA);
            Geometry g2 = reader.read(wktB);

            if (g1 == null || g2 == null || g1.isEmpty() || g2.isEmpty()) {
                return false;
            }
            String t1 = g1.getGeometryType().toUpperCase();
            String t2 = g2.getGeometryType().toUpperCase();
            if (!("POLYGON".equals(t1) || "MULTIPOLYGON".equals(t1))) {
                return false;
            }
            if (!("POLYGON".equals(t2) || "MULTIPOLYGON".equals(t2))) {
                return false;
            }

            // 坐标/几何修复（自相交等）
            if (!hasValidCoordinates(g1)) {
                g1 = g1.buffer(0);
                if (g1.isEmpty())
                    return false;
            }
            if (!hasValidCoordinates(g2)) {
                g2 = g2.buffer(0);
                if (g2.isEmpty())
                    return false;
            }

            // 先做包络矩形快速判断
            Envelope e1 = g1.getEnvelopeInternal();
            Envelope e2 = g2.getEnvelopeInternal();
            if (!e1.intersects(e2)) {
                return false;
            }

            // 使用 PreparedGeometry 提升判断性能
            PreparedGeometry prep = PreparedGeometryFactory.prepare(g1);
            return prep.intersects(g2);
        } catch (ParseException ex) {
            log.error("解析 WKT 失败: {}", ex.getMessage());
            return false;
        } catch (Exception ex) {
            log.error("相交判断出错: {}", ex.getMessage(), ex);
            return false;
        }
    }

    /**
     * 判断两个 WKT（WGS84）是否相等（拓扑）
     * 解析为 `POLYGON`/`MULTIPOLYGON` 后使用 JTS 拓扑相等判断（`equalsTopo`）。
     *
     * @param wktA WKT（WGS84），类型必须为 `POLYGON` 或 `MULTIPOLYGON`
     * @param wktB WKT（WGS84），类型必须为 `POLYGON` 或 `MULTIPOLYGON`
     * @return 是否拓扑相等；解析失败、类型不支持或为空几何返回 false
     *
     * @implNote 坐标非法时先 `buffer(0)` 修复；相等判断不做包络裁剪。
     */
    public boolean equalsWkt(String wktA, String wktB) {
        if (wktA == null || wktA.trim().isEmpty() || wktB == null || wktB.trim().isEmpty()) {
            return false;
        }
        try {
            WKTReader reader = new WKTReader(config.geometryFactory);
            Geometry g1 = reader.read(wktA);
            Geometry g2 = reader.read(wktB);
            if (g1 == null || g2 == null || g1.isEmpty() || g2.isEmpty()) {
                return false;
            }
            String t1 = g1.getGeometryType().toUpperCase();
            String t2 = g2.getGeometryType().toUpperCase();
            if (!("POLYGON".equals(t1) || "MULTIPOLYGON".equals(t1)))
                return false;
            if (!("POLYGON".equals(t2) || "MULTIPOLYGON".equals(t2)))
                return false;
            if (!hasValidCoordinates(g1)) {
                g1 = g1.buffer(0);
                if (g1.isEmpty())
                    return false;
            }
            if (!hasValidCoordinates(g2)) {
                g2 = g2.buffer(0);
                if (g2.isEmpty())
                    return false;
            }
            return g1.equalsTopo(g2);
        } catch (ParseException ex) {
            log.error("解析 WKT 失败: {}", ex.getMessage());
            return false;
        } catch (Exception ex) {
            log.error("相等判断出错: {}", ex.getMessage(), ex);
            return false;
        }
    }

    /**
     * 判断两个 WKT（WGS84）是否脱节（disjoint）
     *
     * @param wktA WKT（WGS84），类型必须为 `POLYGON` 或 `MULTIPOLYGON`
     * @param wktB WKT（WGS84），类型必须为 `POLYGON` 或 `MULTIPOLYGON`
     * @return 是否脱节；解析失败、类型不支持或为空几何返回 false
     *
     * @implNote 先用包络矩形快速裁剪；再用 `PreparedGeometry#disjoint` 判断。
     */
    public boolean disjoint(String wktA, String wktB) {
        if (wktA == null || wktA.trim().isEmpty() || wktB == null || wktB.trim().isEmpty()) {
            return false;
        }
        try {
            WKTReader reader = new WKTReader(config.geometryFactory);
            Geometry g1 = reader.read(wktA);
            Geometry g2 = reader.read(wktB);
            if (g1 == null || g2 == null || g1.isEmpty() || g2.isEmpty())
                return false;
            String t1 = g1.getGeometryType().toUpperCase();
            String t2 = g2.getGeometryType().toUpperCase();
            if (!("POLYGON".equals(t1) || "MULTIPOLYGON".equals(t1)))
                return false;
            if (!("POLYGON".equals(t2) || "MULTIPOLYGON".equals(t2)))
                return false;
            if (!hasValidCoordinates(g1)) {
                g1 = g1.buffer(0);
                if (g1.isEmpty())
                    return false;
            }
            if (!hasValidCoordinates(g2)) {
                g2 = g2.buffer(0);
                if (g2.isEmpty())
                    return false;
            }
            Envelope e1 = g1.getEnvelopeInternal();
            Envelope e2 = g2.getEnvelopeInternal();
            if (!e1.intersects(e2))
                return true;
            PreparedGeometry prep = PreparedGeometryFactory.prepare(g1);
            return prep.disjoint(g2);
        } catch (ParseException ex) {
            log.error("解析 WKT 失败: {}", ex.getMessage());
            return false;
        } catch (Exception ex) {
            log.error("disjoint 判断出错: {}", ex.getMessage(), ex);
            return false;
        }
    }

    /**
     * 判断两个 WKT（WGS84）是否接触（touches）
     * 边界接触但内部不相交。
     *
     * @param wktA WKT（WGS84），类型必须为 `POLYGON` 或 `MULTIPOLYGON`
     * @param wktB WKT（WGS84），类型必须为 `POLYGON` 或 `MULTIPOLYGON`
     * @return 是否接触；解析失败、类型不支持或为空几何返回 false
     */
    public boolean touches(String wktA, String wktB) {
        if (wktA == null || wktA.trim().isEmpty() || wktB == null || wktB.trim().isEmpty())
            return false;
        try {
            WKTReader reader = new WKTReader(config.geometryFactory);
            Geometry g1 = reader.read(wktA);
            Geometry g2 = reader.read(wktB);
            if (g1 == null || g2 == null || g1.isEmpty() || g2.isEmpty())
                return false;
            String t1 = g1.getGeometryType().toUpperCase();
            String t2 = g2.getGeometryType().toUpperCase();
            if (!("POLYGON".equals(t1) || "MULTIPOLYGON".equals(t1)))
                return false;
            if (!("POLYGON".equals(t2) || "MULTIPOLYGON".equals(t2)))
                return false;
            if (!hasValidCoordinates(g1)) {
                g1 = g1.buffer(0);
                if (g1.isEmpty())
                    return false;
            }
            if (!hasValidCoordinates(g2)) {
                g2 = g2.buffer(0);
                if (g2.isEmpty())
                    return false;
            }
            PreparedGeometry prep = PreparedGeometryFactory.prepare(g1);
            return prep.touches(g2);
        } catch (ParseException ex) {
            log.error("解析 WKT 失败: {}", ex.getMessage());
            return false;
        } catch (Exception ex) {
            log.error("touches 判断出错: {}", ex.getMessage(), ex);
            return false;
        }
    }

    /**
     * 判断两个 WKT（WGS84）是否交叉（crosses）
     * 交叉表示几何在维度上“穿过”彼此（多边形相交且交集不是包含/覆盖的关系）。
     *
     * @param wktA WKT（WGS84），类型必须为 `POLYGON` 或 `MULTIPOLYGON`
     * @param wktB WKT（WGS84），类型必须为 `POLYGON` 或 `MULTIPOLYGON`
     * @return 是否交叉；解析失败、类型不支持或为空几何返回 false
     */
    public boolean crosses(String wktA, String wktB) {
        if (wktA == null || wktA.trim().isEmpty() || wktB == null || wktB.trim().isEmpty())
            return false;
        try {
            WKTReader reader = new WKTReader(config.geometryFactory);
            Geometry g1 = reader.read(wktA);
            Geometry g2 = reader.read(wktB);
            if (g1 == null || g2 == null || g1.isEmpty() || g2.isEmpty())
                return false;
            String t1 = g1.getGeometryType().toUpperCase();
            String t2 = g2.getGeometryType().toUpperCase();
            if (!("POLYGON".equals(t1) || "MULTIPOLYGON".equals(t1)))
                return false;
            if (!("POLYGON".equals(t2) || "MULTIPOLYGON".equals(t2)))
                return false;
            if (!hasValidCoordinates(g1)) {
                g1 = g1.buffer(0);
                if (g1.isEmpty())
                    return false;
            }
            if (!hasValidCoordinates(g2)) {
                g2 = g2.buffer(0);
                if (g2.isEmpty())
                    return false;
            }
            PreparedGeometry prep = PreparedGeometryFactory.prepare(g1);
            return prep.crosses(g2);
        } catch (ParseException ex) {
            log.error("解析 WKT 失败: {}", ex.getMessage());
            return false;
        } catch (Exception ex) {
            log.error("crosses 判断出错: {}", ex.getMessage(), ex);
            return false;
        }
    }

    /**
     * 判断 A 是否在 B 内（within）
     *
     * @param wktA WKT（WGS84），A，`POLYGON`/`MULTIPOLYGON`
     * @param wktB WKT（WGS84），B，`POLYGON`/`MULTIPOLYGON`
     * @return A 是否在 B 内；解析失败、类型不支持或为空几何返回 false
     */
    public boolean within(String wktA, String wktB) {
        if (wktA == null || wktA.trim().isEmpty() || wktB == null || wktB.trim().isEmpty())
            return false;
        try {
            WKTReader reader = new WKTReader(config.geometryFactory);
            Geometry g1 = reader.read(wktA);
            Geometry g2 = reader.read(wktB);
            if (g1 == null || g2 == null || g1.isEmpty() || g2.isEmpty())
                return false;
            String t1 = g1.getGeometryType().toUpperCase();
            String t2 = g2.getGeometryType().toUpperCase();
            if (!("POLYGON".equals(t1) || "MULTIPOLYGON".equals(t1)))
                return false;
            if (!("POLYGON".equals(t2) || "MULTIPOLYGON".equals(t2)))
                return false;
            if (!hasValidCoordinates(g1)) {
                g1 = g1.buffer(0);
                if (g1.isEmpty())
                    return false;
            }
            if (!hasValidCoordinates(g2)) {
                g2 = g2.buffer(0);
                if (g2.isEmpty())
                    return false;
            }
            return g1.within(g2);
        } catch (ParseException ex) {
            log.error("解析 WKT 失败: {}", ex.getMessage());
            return false;
        } catch (Exception ex) {
            log.error("within 判断出错: {}", ex.getMessage(), ex);
            return false;
        }
    }

    /**
     * 判断 A 是否包含 B（contains）
     *
     * @param wktA WKT（WGS84），A，`POLYGON`/`MULTIPOLYGON`
     * @param wktB WKT（WGS84），B，`POLYGON`/`MULTIPOLYGON`
     * @return A 是否包含 B；解析失败、类型不支持或为空几何返回 false
     */
    public boolean contains(String wktA, String wktB) {
        if (wktA == null || wktA.trim().isEmpty() || wktB == null || wktB.trim().isEmpty())
            return false;
        try {
            WKTReader reader = new WKTReader(config.geometryFactory);
            Geometry g1 = reader.read(wktA);
            Geometry g2 = reader.read(wktB);
            if (g1 == null || g2 == null || g1.isEmpty() || g2.isEmpty())
                return false;
            String t1 = g1.getGeometryType().toUpperCase();
            String t2 = g2.getGeometryType().toUpperCase();
            if (!("POLYGON".equals(t1) || "MULTIPOLYGON".equals(t1)))
                return false;
            if (!("POLYGON".equals(t2) || "MULTIPOLYGON".equals(t2)))
                return false;
            if (!hasValidCoordinates(g1)) {
                g1 = g1.buffer(0);
                if (g1.isEmpty())
                    return false;
            }
            if (!hasValidCoordinates(g2)) {
                g2 = g2.buffer(0);
                if (g2.isEmpty())
                    return false;
            }
            return g1.contains(g2);
        } catch (ParseException ex) {
            log.error("解析 WKT 失败: {}", ex.getMessage());
            return false;
        } catch (Exception ex) {
            log.error("contains 判断出错: {}", ex.getMessage(), ex);
            return false;
        }
    }

    /**
     * 判断两个 WKT（WGS84）是否重叠（overlaps）
     * 多边形维度一致，交集维度同为 2，且交集既不是包含也不是相等。
     *
     * @param wktA WKT（WGS84），类型必须为 `POLYGON` 或 `MULTIPOLYGON`
     * @param wktB WKT（WGS84），类型必须为 `POLYGON` 或 `MULTIPOLYGON`
     * @return 是否重叠；解析失败、类型不支持或为空几何返回 false
     */
    public boolean overlaps(String wktA, String wktB) {
        if (wktA == null || wktA.trim().isEmpty() || wktB == null || wktB.trim().isEmpty())
            return false;
        try {
            WKTReader reader = new WKTReader(config.geometryFactory);
            Geometry g1 = reader.read(wktA);
            Geometry g2 = reader.read(wktB);
            if (g1 == null || g2 == null || g1.isEmpty() || g2.isEmpty())
                return false;
            String t1 = g1.getGeometryType().toUpperCase();
            String t2 = g2.getGeometryType().toUpperCase();
            if (!("POLYGON".equals(t1) || "MULTIPOLYGON".equals(t1)))
                return false;
            if (!("POLYGON".equals(t2) || "MULTIPOLYGON".equals(t2)))
                return false;
            if (!hasValidCoordinates(g1)) {
                g1 = g1.buffer(0);
                if (g1.isEmpty())
                    return false;
            }
            if (!hasValidCoordinates(g2)) {
                g2 = g2.buffer(0);
                if (g2.isEmpty())
                    return false;
            }
            PreparedGeometry prep = PreparedGeometryFactory.prepare(g1);
            return prep.overlaps(g2);
        } catch (ParseException ex) {
            log.error("解析 WKT 失败: {}", ex.getMessage());
            return false;
        } catch (Exception ex) {
            log.error("overlaps 判断出错: {}", ex.getMessage(), ex);
            return false;
        }
    }

    /**
     * 生成轨迹轮廓（不拆分、不裁剪、不平滑）
     * 根据轨迹点与总宽度，使用线缓冲构建单个轮廓（Polygon），并返回包含亩数、WKT、时间范围与点集的 `OutlinePart`。
     * 步骤：过滤/排序轨迹点 → 在米制投影下线简化+线缓冲 → 若为MultiPolygon仅做轻微开运算保留缝隙并选面积最大面 → 计算亩数与 WKT。
     *
     * @param seg         轨迹点列表（WGS84），至少 3 个有效点
     * @param totalWidthM 总宽度（米，左右合计），必须非负
     * @return `OutlinePart`：结果几何为 `Polygon`，包含 `mu`、`wkt`、起止时间与有效点集
     * @throws IllegalArgumentException 如果参数无效或轨迹点不足
     * @throws Exception                如果几何运算或坐标转换失败
     */
    public OutlinePart getOutline(List<TrackPoint> seg, double totalWidthM) throws Exception {
        long t0 = System.currentTimeMillis();

        // 参数校验
        if (seg == null) {
            throw new IllegalArgumentException("轨迹点列表不能为空");
        }
        if (totalWidthM < 0) {
            throw new IllegalArgumentException("总宽度必须为非负数");
        }

        // 计算半宽用于缓冲
        double widthM = totalWidthM / 2.0;

        // 步骤1: 过滤异常点（越界与(0,0)），按时间升序排序
        List<TrackPoint> points = filterAndSortTrackPoints(seg);
        log.trace("[getOutline] 轨迹点过滤完成，原始点数: {}, 过滤后点数: {}", seg.size(), points.size());

        // 校验有效点数
        if (points.size() < 3) {
            throw new IllegalArgumentException("轨迹段至少需要3个有效点");
        }

        // 步骤2: 提取时间范围（忽略空时间）
        LocalDateTime startTime = null;
        LocalDateTime endTime = null;
        for (TrackPoint p : points) {
            if (p.getTime() != null) {
                if (startTime == null) {
                    startTime = p.getTime();
                }
                endTime = p.getTime(); // 总是更新为最后一个有效时间
            }
        }

        // 步骤3: 构建线缓冲轮廓（WGS84坐标系）
        long tBuildStart = System.currentTimeMillis();
        Geometry outline = buildSimpleLineBuffer(points, widthM);
        long tBuildEnd = System.currentTimeMillis();
        log.trace("[getOutline] 轮廓构建完成（线缓冲），耗时: {}ms", (tBuildEnd - tBuildStart));

        // 步骤4: 计算几何属性
        double mu = calcMu(outline); // 计算面积（亩）
        String wkt = toWkt(outline); // 转换为WKT格式

        // 记录最终结果信息
        log.debug("[getOutline] 返回结果 type={} 面积={}亩 点数={}",
                outline.getGeometryType(), mu, points.size());

        long t1 = System.currentTimeMillis();
        log.debug("[getOutline] 总耗时={}ms", (t1 - t0));

        // 构建并返回结果对象
        OutlinePart op = new OutlinePart();
        op.setOutline(outline);
        op.setStartTime(startTime);
        op.setEndTime(endTime);
        op.setMu(mu);
        op.setWkt(wkt);
        op.setTrackPoints(points);
        op.setTotalWidthM(totalWidthM);
        return op;
    }

    public SplitRoadResult splitRoad(List<TrackPoint> wgs84Points, double totalWidthM) {
        return splitRoad(wgs84Points, totalWidthM, config.DEFAULT_MAX_OUTLINE_SEGMENTS);
    }

    public SplitRoadResult splitRoad(List<TrackPoint> wgs84Points, double totalWidthM, int maxSegments) {
        SplitRoadResult result = new SplitRoadResult();
        result.setTotalWidthM(totalWidthM);
        result.setOutline(config.EMPTY_GEOMETRY);
        result.setWkt(config.EMPTY_WKT);
        result.setMu(0);

        // 参数校验，seg必须大于3个点
        if (CollUtil.isEmpty(wgs84Points) || wgs84Points.size() < 3) {
            log.warn("轨迹点列表必须包含至少3个有效点");
            return result;
        }
        // 参数校验，totalWidthM必须大于1米
        if (totalWidthM < 1) {
            log.warn("作业总宽度不能小于1米");
            return result;
        }
        // 参数校验，maxSegments必须大于0
        if (maxSegments <= 0) {
            log.warn("返回的区块数量上限必须大于0");
            return result;
        }

        log.debug("参数 wgs84Points.size={}, totalWidthM={}, maxSegments={}", wgs84Points.size(), totalWidthM,
                maxSegments);
        // 循环seg，如果时间为null，抛弃该点
        wgs84Points.removeIf(p -> p.getTime() == null);
        // 先将seg进行按时间升序排序
        wgs84Points.sort(Comparator.comparing(TrackPoint::getTime));
        // 过滤经纬度异常的点（纬度[-90,90]，经度[-180,180]）
        wgs84Points.removeIf(p -> p.getLat() < -90 || p.getLat() > 90 || p.getLon() < -180 || p.getLon() > 180);
        // 过滤经纬度为0的点
        wgs84Points.removeIf(p -> p.getLat() == 0 && p.getLon() == 0);
        if (CollUtil.isEmpty(wgs84Points) || wgs84Points.size() < 3) {
            log.warn("过滤后轨迹点列表不足3个点");
            return result;
        }
        log.debug("过滤后 wgs84Points.size={}", wgs84Points.size());

        // 将wgs84Points进行wgs84到高斯投影的转换
        List<TrackPoint> gaussPoints = convertWgs84ToGaussProjection(wgs84Points);
        if (CollUtil.isEmpty(gaussPoints) || gaussPoints.size() < 3) {
            log.warn("转换后的轨迹点列表不足3个有效点");
            return result;
        }

        // 过滤速度范围的点（高斯投影坐标系，米单位）
        /*
         * gaussPoints = filterBySpeedRangeMeter(gaussPoints, config.MIN_WORK_SPEED_KMH,
         * config.WORK_MAX_SPEED_KMH);
         * if (CollUtil.isEmpty(gaussPoints) || gaussPoints.size() < 3) {
         * log.warn("速度过滤后轨迹点列表不足3个有效点");
         * return result;
         * }
         */

        // 计算单侧缓冲宽度（totalWidthM是总宽度，需要除以2）
        double halfWidth = totalWidthM / 2.0;
        double halfHalfWidth = halfWidth / 2.0;

        // 实现按上报频率分组：统计速度过滤后的轨迹点能分出多少频率相同的组
        // gaussPoints已经在前面按时间排序（第2658行），无需再次排序
        if (gaussPoints.size() < 2) {
            log.warn("速度过滤后轨迹点数不足2个，无法进行频率分组");
            return result;
        }

        log.info("速度过滤后轨迹点数: {}", gaussPoints.size());

        // 第一步：定义标准频率模式（1秒、5秒、10秒、15秒）
        int[] standardFrequencies = { 1, 5, 10, 15 };
        Map<Integer, Integer> timeIntervalStats = new LinkedHashMap<>(); // 标准频率 -> 出现次数

        for (int i = 1; i < gaussPoints.size(); i++) {
            long timeDiff = java.time.Duration.between(gaussPoints.get(i - 1).getTime(), gaussPoints.get(i).getTime())
                    .getSeconds();

            // 只处理1-17秒的间隔（15秒+2秒容错）
            if (timeDiff >= 1 && timeDiff <= 17) {
                int standardFreq = mapToStandardFrequency((int) timeDiff, standardFrequencies);
                if (standardFreq != -1) {
                    timeIntervalStats.merge(standardFreq, 1, Integer::sum);
                }
            }
        }

        if (timeIntervalStats.isEmpty()) {
            log.warn("未找到符合标准频率的时间间隔");
            return result;
        }

        // 按标准频率排序输出
        List<Integer> sortedFrequencies = new ArrayList<>(timeIntervalStats.keySet());
        Collections.sort(sortedFrequencies);

        log.info("发现的标准频率模式（共{}种）:", sortedFrequencies.size());
        int totalIntervals = timeIntervalStats.values().stream().mapToInt(Integer::intValue).sum();

        for (int freq : sortedFrequencies) {
            int count = timeIntervalStats.get(freq);
            double percentage = (count * 100.0) / totalIntervals;
            log.info("频率{}秒: 出现{}次 ({}%)", freq, count, percentage);
        }

        // 第二步：根据发现的频率模式进行分组
        Map<Integer, List<List<TrackPoint>>> frequencyGroups = new LinkedHashMap<>();
        List<TrackPoint> currentGroup = new ArrayList<>();
        int currentFrequency = -1;

        for (int i = 0; i < gaussPoints.size(); i++) {
            TrackPoint currentPoint = gaussPoints.get(i);

            if (i == 0) {
                currentGroup.add(currentPoint);
                continue;
            }

            // 计算与前一个点的时间间隔（秒）
            long timeDiff = Duration.between(gaussPoints.get(i - 1).getTime(), currentPoint.getTime())
                    .getSeconds();

            // 如果时间间隔超过15秒，作为异常间隔处理（结束当前组）
            if (timeDiff > 15) {
                log.warn("发现时间间隔{}秒超过15秒，作为异常间隔处理", timeDiff);
                // 异常间隔，结束当前组
                if (currentGroup.size() >= 2) {
                    frequencyGroups.computeIfAbsent(currentFrequency, k -> new ArrayList<>())
                            .add(new ArrayList<>(currentGroup));
                }
                currentGroup.clear();
                currentGroup.add(currentPoint);
                currentFrequency = -1;
                continue;
            }

            // 将实际间隔映射到标准频率（正负2秒容错）
            int standardFreq = mapToStandardFrequency((int) timeDiff, standardFrequencies);

            if (standardFreq != -1) {
                if (currentFrequency == -1) {
                    currentFrequency = standardFreq;
                }

                if (standardFreq == currentFrequency) {
                    currentGroup.add(currentPoint);
                } else {
                    // 频率变化，结束当前组
                    if (currentGroup.size() >= 2) {
                        frequencyGroups.computeIfAbsent(currentFrequency, k -> new ArrayList<>())
                                .add(new ArrayList<>(currentGroup));
                    }
                    currentGroup.clear();
                    currentGroup.add(currentPoint);
                    currentFrequency = standardFreq;
                }
            } else {
                // 未识别的间隔（超过15秒或不在标准频率范围内），结束当前组并丢弃这组数据
                if (currentGroup.size() >= 2) {
                    frequencyGroups.computeIfAbsent(currentFrequency, k -> new ArrayList<>())
                            .add(new ArrayList<>(currentGroup));
                }
                currentGroup.clear();
                currentGroup.add(currentPoint);
                currentFrequency = -1;
            }
        }

        // 处理最后一组
        if (currentGroup.size() >= 2 && currentFrequency != -1) {
            frequencyGroups.computeIfAbsent(currentFrequency, k -> new ArrayList<>())
                    .add(new ArrayList<>(currentGroup));
        }

        int totalGroups = 0;

        // 按标准频率顺序输出
        for (int freq : standardFrequencies) {
            if (!frequencyGroups.containsKey(freq)) {
                continue; // 跳过没有数据的频率
            }

            List<List<TrackPoint>> groups = frequencyGroups.get(freq);
            int groupCount = groups.size();
            int pointCount = groups.stream().mapToInt(List::size).sum();

            totalGroups += groupCount;

            log.info("频率{}秒: {}个组, {}个点", freq, groupCount, pointCount);

            // 打印每个组的详细信息
            for (int i = 0; i < groups.size(); i++) {
                List<TrackPoint> group = groups.get(i);
                if (!group.isEmpty()) {
                    log.debug("  组{}: {}个点, 时间:{} ~ {}",
                            i + 1,
                            group.size(),
                            group.get(0).getTime(),
                            group.get(group.size() - 1).getTime());
                }
            }
        }

        log.info("总频率模式数: {}", frequencyGroups.size());
        log.info("总组数: {}", totalGroups);

        // 获取最小间隔的分组
        List<Geometry> allPolygons = new ArrayList<>();
        if (!frequencyGroups.isEmpty()) {
            int minFrequency = frequencyGroups.keySet().stream().min(Integer::compareTo).orElse(0);
            List<List<TrackPoint>> minFreqGroups = frequencyGroups.get(minFrequency);
            int minFreqGroupCount = minFreqGroups.size();
            int minFreqPointCount = minFreqGroups.stream().mapToInt(List::size).sum();
            log.info("最小间隔{}秒: {}个组, {}个点", minFrequency, minFreqGroupCount, minFreqPointCount);

            for (int k = 0; k < minFreqGroups.size(); k++) {
                List<TrackPoint> group = minFreqGroups.get(k);
                // 计算每一组的距离中位数
                if (group.size() >= 3) {
                    List<Double> distances = new ArrayList<>();
                    for (int j = 1; j < group.size(); j++) {
                        TrackPoint prev = group.get(j - 1);
                        TrackPoint curr = group.get(j);
                        double distance = Math.hypot(
                                curr.getLon() - prev.getLon(),
                                curr.getLat() - prev.getLat());
                        distances.add(distance);
                    }

                    if (!distances.isEmpty()) {
                        distances.sort(Double::compareTo);
                        double medianDistance;
                        if (distances.size() % 2 == 0) {
                            medianDistance = (distances.get(distances.size() / 2 - 1) +
                                    distances.get(distances.size() / 2)) / 2.0;
                        } else {
                            medianDistance = distances.get(distances.size() / 2);
                        }

                        // 计算平均数、最大距离、最小距离
                        double sumDistance = distances.stream().mapToDouble(Double::doubleValue).sum();
                        double avgDistance = sumDistance / distances.size();
                        double maxDistance = distances.get(distances.size() - 1);

                        log.info("  组{}: {}个点, 距离统计: 中位数={}米, 平均数={}米, 最大={}米",
                                k + 1, group.size(), medianDistance, avgDistance, maxDistance);

                        if (0 < medianDistance && medianDistance <= 20) {
                            try {
                                double segmentDistanceThreshold = Math.max(
                                        config.MIN_SEGMENT_DISTANCE_THRESHOLD_MAP.floorEntry(totalWidthM).getValue(),
                                        minFrequency * totalWidthM);

                                List<List<TrackPoint>> segments = new ArrayList<>();
                                List<TrackPoint> currentSegment = new ArrayList<>();

                                log.debug("开始按{}米距离阈值对组{}进行轨迹分段", segmentDistanceThreshold, k + 1);

                                // 使用当前组的点进行分段，而不是所有轨迹点
                                for (int i = 0; i < group.size(); i++) {
                                    TrackPoint currentPoint = group.get(i);

                                    if (currentSegment.isEmpty()) {
                                        // 第一段或新段开始的第一个点
                                        currentSegment.add(currentPoint);
                                        continue;
                                    }

                                    // 计算当前点与段内最后一个点的距离
                                    TrackPoint lastPoint = currentSegment.get(currentSegment.size() - 1);
                                    double distance = Math.hypot(
                                            currentPoint.getLon() - lastPoint.getLon(),
                                            currentPoint.getLat() - lastPoint.getLat());

                                    if (distance > segmentDistanceThreshold) {
                                        // 距离超过阈值，结束当前段，开始新段
                                        if (currentSegment.size() >= 2) {
                                            segments.add(new ArrayList<>(currentSegment));
                                            log.debug("组{}分段完成: 段内点数={}, 最后距离={}m", k + 1, currentSegment.size(),
                                                    distance);
                                        }
                                        currentSegment.clear();
                                    }

                                    currentSegment.add(currentPoint);
                                }

                                // 处理最后一段
                                if (currentSegment.size() >= 2) {
                                    segments.add(currentSegment);
                                    log.debug("组{}最后一段: 段内点数={}", k + 1, currentSegment.size());
                                }

                                log.debug("组{}轨迹分段完成: 共{}段", k + 1, segments.size());

                                // 对每一段分别构建线缓冲
                                for (int segIndex = 0; segIndex < segments.size(); segIndex++) {
                                    List<TrackPoint> segmentPoints = segments.get(segIndex);

                                    if (segmentPoints.size() < 2) {
                                        log.warn("组{}第{}段点数不足2个，跳过", k + 1, segIndex);
                                        continue;
                                    }

                                    // 构建当前段的线几何
                                    Coordinate[] coords = new Coordinate[segmentPoints.size()];
                                    for (int i = 0; i < segmentPoints.size(); i++) {
                                        TrackPoint pt = segmentPoints.get(i);
                                        coords[i] = new Coordinate(pt.getLon(), pt.getLat());
                                    }

                                    log.debug("有 {} 个点需要构建线缓冲", coords.length);

                                    // 对大型轨迹段进行分组处理（每500个点一组）
                                    Geometry buffer;
                                    if (coords.length > config.BATCH_SIZE) {
                                        log.debug("组{}第{}段点数过多，进行分组处理", k + 1, segIndex);

                                        // 计算需要多少个批次
                                        final int totalCoords = coords.length;
                                        final int overlapSize = 1;
                                        final int batchSize = config.BATCH_SIZE;
                                        final GeometryFactory geometryFactory = config.geometryFactory;
                                        final double hw = halfWidth;
                                        final int segmentIndex = segIndex;
                                        final int groupIndex = k + 1;

                                        // 使用Java 8 Stream API并行处理批次
                                        List<Geometry> batchBuffers = IntStream.range(0, totalCoords)
                                                .filter(i -> i % (batchSize - overlapSize) == 0)
                                                .parallel() // 并行执行以充分利用CPU核心
                                                .mapToObj(i -> {
                                                    int endIndex = Math.min(i + batchSize, totalCoords);
                                                    Coordinate[] batchCoords = new Coordinate[endIndex - i];
                                                    System.arraycopy(coords, i, batchCoords, 0, endIndex - i);

                                                    try {
                                                        LineString batchLine = geometryFactory
                                                                .createLineString(batchCoords);
                                                        Geometry batchBuffer = batchLine.buffer(hw);

                                                        log.debug("组{}第{}段批次{}处理完成: 点数={}, 缓冲类型={}",
                                                                groupIndex, segmentIndex,
                                                                i / (batchSize - overlapSize) + 1,
                                                                batchCoords.length, batchBuffer.getGeometryType());

                                                        return !batchBuffer.isEmpty() ? batchBuffer : null;
                                                    } catch (Exception e) {
                                                        log.warn("批次处理失败 (i={}, endIndex={}): {}", i, endIndex,
                                                                e.getMessage());
                                                        return null;
                                                    }
                                                })
                                                .filter(geo -> geo != null)
                                                .collect(Collectors.toList());

                                        // 在缓冲结果处理前进行合并操作，模仿后续的合并逻辑
                                        if (batchBuffers.isEmpty()) {
                                            buffer = config.EMPTY_GEOMETRY;
                                        } else {
                                            try {
                                                // 先对每个多边形应用缓冲扩展，增加合并的敏感度
                                                List<Geometry> bufferedPolygons = batchBuffers.parallelStream() // 并行处理膨胀操作
                                                        .map(polygon -> {
                                                            try {
                                                                return polygon.buffer(halfHalfWidth);
                                                            } catch (Exception e) {
                                                                log.warn("多边形膨胀处理失败: {}", e.getMessage());
                                                                return polygon; // 失败时返回原始多边形
                                                            }
                                                        })
                                                        .filter(geo -> geo != null)
                                                        .collect(Collectors.toList());

                                                UnaryUnionOp unionOp = new UnaryUnionOp(bufferedPolygons);
                                                buffer = unionOp.union();

                                                // 合并后再收缩回原始大小（减去缓冲距离）
                                                if (!buffer.isEmpty()) {
                                                    buffer = buffer.buffer(-halfHalfWidth);
                                                }
                                            } catch (Exception e) {
                                                log.warn("批次缓冲合并失败: {}", e.getMessage());
                                                // 如果合并失败，尝试返回第一个有效几何
                                                buffer = batchBuffers.get(0);
                                            }
                                        }
                                    } else {
                                        // 点数较少时直接处理
                                        LineString trackLine = config.geometryFactory.createLineString(coords);
                                        buffer = trackLine.buffer(halfWidth);
                                    }

                                    log.debug("组{}第{}段线缓冲构建完成: 原始点数={}, 缓冲类型={}, 顶点数={}",
                                            k + 1, segIndex, segmentPoints.size(), buffer.getGeometryType(),
                                            buffer.getNumPoints());

                                    // 处理缓冲结果
                                    if (buffer.isEmpty()) {
                                        log.warn("组{}第{}段线缓冲为空", k + 1, segIndex);
                                    } else if (buffer instanceof Polygon) {
                                        // 单多边形
                                        allPolygons.add(buffer);
                                        log.debug("组{}第{}段添加单多边形", k + 1, segIndex);
                                    } else if (buffer instanceof MultiPolygon) {
                                        // 多多边形
                                        MultiPolygon multiPoly = (MultiPolygon) buffer;
                                        for (int i = 0; i < multiPoly.getNumGeometries(); i++) {
                                            allPolygons.add(multiPoly.getGeometryN(i));
                                            log.debug("组{}第{}段添加子多边形: 索引={}", k + 1, segIndex, i);
                                        }
                                        log.debug("组{}第{}段添加多多边形: 子多边形数量={}", k + 1, segIndex,
                                                multiPoly.getNumGeometries());
                                    } else {
                                        log.warn("组{}第{}段未知的缓冲几何类型: {}", k + 1, segIndex, buffer.getGeometryType());
                                    }
                                }

                                log.debug("组{}所有段处理完成: 总段数={}, 生成子多边形数量={}", k + 1, segments.size(),
                                        allPolygons.size());
                            } catch (Exception e) {
                                log.warn("组{}线缓冲构建失败: {}", k + 1, e.getMessage());
                            }
                        }
                    }
                }
            }
        }

        // 将allPolygons再次进行一次合并
        List<OutlinePart> parts = new ArrayList<>();
        if (CollUtil.isNotEmpty(allPolygons)) {
            try {
                // 使用UnaryUnionOp合并相邻或相交的多边形
                // 先对每个多边形应用缓冲扩展，增加合并的敏感度
                List<Geometry> bufferedPolygons = new ArrayList<>();
                for (Geometry polygon : allPolygons) {
                    // 使用配置中的合并缓冲距离扩展多边形边界
                    Geometry bufferedPolygon = polygon.buffer(halfHalfWidth);
                    bufferedPolygons.add(bufferedPolygon);
                }

                UnaryUnionOp unionOp = new UnaryUnionOp(bufferedPolygons);
                Geometry mergedGeometry = unionOp.union();

                // 合并后再收缩回原始大小（减去缓冲距离）
                if (!mergedGeometry.isEmpty()) {
                    mergedGeometry = mergedGeometry.buffer(-halfHalfWidth);
                }

                // 将合并后的几何转换为几何列表（避免类型强转）
                List<Geometry> mergedGeometries = new ArrayList<>();
                if (mergedGeometry instanceof Polygon || mergedGeometry instanceof MultiPolygon) {
                    if (mergedGeometry instanceof Polygon) {
                        mergedGeometries.add(mergedGeometry);
                    } else {
                        MultiPolygon multiPoly = (MultiPolygon) mergedGeometry;
                        for (int i = 0; i < multiPoly.getNumGeometries(); i++) {
                            mergedGeometries.add(multiPoly.getGeometryN(i));
                        }
                    }
                }

                // 为每个合并后的几何创建OutlinePart
                for (Geometry geom : mergedGeometries) {
                    OutlinePart part = new OutlinePart();
                    part.setTotalWidthM(totalWidthM);
                    part.setOutline(geom);
                    part.setWkt(gaussGeometryToWgs84Wkt(geom, wgs84Points.get(0).getLon()));
                    // 用当前段的所有轨迹点过滤
                    List<TrackPoint> geometryGaussPoints = filterGaussPointsByGeometry(gaussPoints, geom);
                    part.setTrackPoints(
                            gaussPointsToWgs84(geometryGaussPoints, wgs84Points.get(0).getLon()));
                    if (!geometryGaussPoints.isEmpty()) {
                        part.setStartTime(
                                geometryGaussPoints.stream().map(TrackPoint::getTime).min(LocalDateTime::compareTo)
                                        .orElse(null));
                        part.setEndTime(
                                geometryGaussPoints.stream().map(TrackPoint::getTime).max(LocalDateTime::compareTo)
                                        .orElse(null));
                    }
                    part.setMu(calcMuByWgs84Wkt(part.getWkt()));
                    parts.add(part);
                }

                log.info("多边形合并完成: 原始多边形数量={}, 合并后几何数量={}",
                        allPolygons.size(), mergedGeometries.size());

            } catch (Exception e) {
                log.warn("多边形合并失败，使用原始多边形: {}", e.getMessage());
                // 如果合并失败，回退到原始逻辑
                for (Geometry geom : allPolygons) {
                    OutlinePart part = new OutlinePart();
                    part.setTotalWidthM(totalWidthM);
                    part.setOutline(geom);
                    part.setWkt(gaussGeometryToWgs84Wkt(geom, wgs84Points.get(0).getLon()));
                    // 用当前段的所有轨迹点过滤
                    List<TrackPoint> geometryGaussPoints = filterGaussPointsByGeometry(gaussPoints, geom);
                    part.setTrackPoints(
                            gaussPointsToWgs84(geometryGaussPoints, wgs84Points.get(0).getLon()));
                    if (!geometryGaussPoints.isEmpty()) {
                        part.setStartTime(
                                geometryGaussPoints.stream().map(TrackPoint::getTime).min(LocalDateTime::compareTo)
                                        .orElse(null));
                        part.setEndTime(
                                geometryGaussPoints.stream().map(TrackPoint::getTime).max(LocalDateTime::compareTo)
                                        .orElse(null));
                    }
                    part.setMu(calcMuByWgs84Wkt(part.getWkt()));
                    parts.add(part);
                }
            }
        }

        if (CollUtil.isNotEmpty(parts)) {
            double mu = 0.0;

            // 构建总的MultiPolygon outline
            List<Polygon> polygons = new ArrayList<>();
            for (OutlinePart op : parts) {
                mu += op.getMu();
                polygons.add((Polygon) op.getOutline());
            }

            // 创建MultiPolygon
            Geometry finalOutline;
            if (polygons.size() == 1) {
                finalOutline = polygons.get(0);
            } else {
                finalOutline = config.geometryFactory.createMultiPolygon(polygons.toArray(new Polygon[0]));
            }

            result.setParts(parts);
            result.setOutline(finalOutline);
            result.setWkt(gaussGeometryToWgs84Wkt(finalOutline, wgs84Points.get(0).getLon()));
            result.setMu(mu);
            // 使用parts列表中所有OutlinePart的最早和最晚时间
            result.setStartTime(
                    parts.stream()
                            .map(OutlinePart::getStartTime)
                            .filter(startTime -> startTime != null)
                            .min(LocalDateTime::compareTo)
                            .orElse(null));
            result.setEndTime(
                    parts.stream()
                            .map(OutlinePart::getEndTime)
                            .filter(endTime -> endTime != null)
                            .max(LocalDateTime::compareTo)
                            .orElse(null));
        }

        return result;
    }

}