package sunyu.util;

import java.time.LocalDateTime;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.concurrent.ConcurrentHashMap;
import java.util.stream.Collectors;

import org.geotools.geometry.jts.JTS;
import org.geotools.referencing.CRS;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Envelope;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.LineString;
import org.locationtech.jts.geom.MultiLineString;
import org.locationtech.jts.geom.MultiPolygon;
import org.locationtech.jts.geom.Polygon;
import org.locationtech.jts.geom.prep.PreparedGeometry;
import org.locationtech.jts.geom.prep.PreparedGeometryFactory;
import org.locationtech.jts.io.ParseException;
import org.locationtech.jts.io.WKTReader;
import org.locationtech.jts.operation.buffer.BufferOp;
import org.locationtech.jts.operation.buffer.BufferParameters;
import org.locationtech.jts.operation.union.UnaryUnionOp;
import org.locationtech.jts.simplify.DouglasPeuckerSimplifier;
import org.locationtech.jts.simplify.TopologyPreservingSimplifier;
import org.opengis.referencing.crs.CoordinateReferenceSystem;
import org.opengis.referencing.operation.MathTransform;

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
 * 公共方法一览：
 * - builder()：创建构建器，用于配置并构建 `GisUtil`
 * - close()：释放资源（实现 `AutoCloseable`）
 * - toWkt(Geometry)：将几何转换为 WKT（统一 WGS84），含坐标识别与修复
 * - fromWkt(String)：解析 WKT（WGS84）为 Geometry 并转换到高斯-克吕格米制坐标
 * - haversine(CoordinatePoint, CoordinatePoint)：计算两点大圆距离（米，WGS84）
 * - calcMu(Geometry)：计算几何面积的亩数（球面公式）
 * - calcMu(String)：解析 WKT 后计算亩数
 * - intersection(String, String)：计算两 WKT 相交的几何与面积，返回 `WktIntersectionResult`
 * - intersects(String, String)：判断两 WKT 是否相交（WGS84，支持 `POLYGON`/`MULTIPOLYGON`）
 * - equalsWkt(String, String)：判断两 WKT 是否拓扑相等
 * - disjoint(String, String)：判断是否脱节
 * - touches(String, String)：判断是否接触（边界接触）
 * - crosses(String, String)：判断是否交叉
 * - within(String, String)：判断 A 是否在 B 内
 * - contains(String, String)：判断 A 是否包含 B
 * - overlaps(String, String)：判断是否重叠
 * - pointInPolygon(CoordinatePoint, String)：判断点是否在多边形内（含边界）
 * - getOutline(List<TrackPoint>, double)：生成轨迹轮廓（`Polygon`），返回 `OutlinePart`
 * - splitRoad(List<TrackPoint>, double)：按总宽度对轨迹进行分段并返回结果
 * - splitRoad(List<TrackPoint>, double, Integer)：指定最大段数的道路分段
 *
 * 概览：
 * - 坐标系：默认 WGS84；轮廓构建与形态学在高斯-克吕格米制投影（按首点经度分带）进行，结果回转 WGS84。
 * - 常量：距离计算使用平均地球半径 R=6371000；面积计算在 ringArea 使用 WGS84 半径 6378137（与 Turf.js 对齐）。
 * - 轮廓：线缓冲构建（折线简化+缓冲，拐角/端头由 BufferParameters 控制）；getOutline 不切割、不平滑；splitRoad
 * 可选轻微开运算、外缘细长裁剪与缝隙蚀刻。
 * - 边界：`pointInPolygon` 边界视为内；几何谓词遵循 JTS 语义。
 *
 * 设计要点：
 * - 坐标系：内部优先在 WGS84 下处理；涉及形态学与面积计算时使用高斯-克吕格米制投影（6 度分带）。
 * - 变换缓存：按分带缓存 CRS 与 MathTransform，避免重复构建。
 * - 轮廓构建：线简化+线缓冲并合并，控制几何规模与性能。
 * - 清理过滤：按面积（亩数）等规则过滤过小碎片（仅 splitRoad 中使用）。
 * - 输出：支持 WKT 输出并在必要时进行坐标系识别与修复。
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
     * 持有常量、默认参数、线程安全缓存与几何工厂；供 `GisUtil` 使用。
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

        // 默认轮廓返回的最多多边形数量（TopN）
        private final int DEFAULT_MAX_OUTLINE_SEGMENTS = 10;

        // 圆近似细分：4，提升性能，误差满足道路宽度场景
        private final int DEFAULT_BUFFER_QUADRANT = 2;

        // 作业最大速度阈值（km/h），用于前置速度过滤；可通过 Builder 配置
        private double WORK_MAX_SPEED_KMH = 20.0;
        // 作业最小速度阈值（km/h），用于前置速度过滤；可通过 Builder 配置
        private double MIN_WORK_SPEED_KMH = 1.0;
        // 最小亩数动态阈值（亩），用于 splitRoad 动态过滤小块
        private double MIN_MU_DYNAMIC_THRESHOLD_MU = 0.5;

        // 外缘细长裁剪开关（仅 splitRoad 使用；getOutline 不使用）
        private boolean ENABLE_OUTER_THIN_TRIM = true;
        // 外缘细长裁剪半径系数（相对单侧宽度）
        private double THIN_TRIM_RADIUS_FACTOR = 1.6;

        // 线段断裂距离系数（倍数*单侧宽度），超过则切分会话
        private double LINE_BREAK_FACTOR = 4.0;
        // 线简化公差系数（倍数*单侧宽度），用于Douglas-Peucker
        private double LINE_SIMPLIFY_TOL_FACTOR = 0.1;

        // 凹包开运算（保留缝隙/洞）：在米制下做轻微开运算
        private boolean ENABLE_CONCAVE_SMOOTHING = true;
        // 开运算半径系数（相对于单侧宽度widthM）；建议 0.3–0.8
        private double CONCAVE_SMOOTH_RADIUS_FACTOR = 0.4;

        // 线缓冲样式（更锐利边界）：拐角样式、端头样式与mitre限制
        private int BUFFER_JOIN_STYLE = BufferParameters.JOIN_MITRE;
        private int BUFFER_END_CAP_STYLE = BufferParameters.CAP_FLAT;
        private double BUFFER_MITRE_LIMIT = 2.0;

        // 可选：轻微蚀刻扩大缝隙（在米制下对并集做负缓冲）
        private boolean ENABLE_GAP_ENHANCE_ERODE = true;
        private double GAP_ENHANCE_ERODE_FACTOR = 0.2;

        // 直线道路段断开控制（防止窄直线连接合并两块大轮廓）
        private boolean ENABLE_ROAD_STRAIGHT_BREAK = true;
        private double ROAD_STRAIGHT_DIFF_THRESHOLD_DEG = 6.0;
        private int ROAD_STRAIGHT_WINDOW_POINTS = 8;
        private double ROAD_BREAK_DIST_FACTOR = 2.0;

    }

    /**
     * 构建器类，用于构建GisUtil实例
     * 实现构建器模式，允许灵活配置和创建GisUtil对象
     */
    public static class Builder {
        // 配置对象，包含各种常量和默认值
        private final Config config = new Config();

        /**
         * 构建GisUtil实例
         * 使用预配置的参数创建GisUtil对象
         *
         * @return GisUtil对象，已初始化完成可直接使用
         */
        public GisUtil build() {
            return new GisUtil(config);
        }

    }

    /**
     * 关闭资源（AutoCloseable接口实现）
     * 用于释放GIS工具类使用的各种资源，如文件句柄、网络连接等
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
     * 获取/缓存坐标参考系（CRS）。
     * 使用并发缓存按 code（如 EPSG:4326）惰性解析并存储，避免重复调用 CRS.decode。
     *
     * @param code 坐标系代码，如 EPSG:4326
     * @return 解析后的CRS，线程安全缓存
     */
    private CoordinateReferenceSystem getCachedCRS(String code) {
        return config.crsCache.computeIfAbsent(code, key -> {
            try {
                return CRS.decode(key, true);
            } catch (Exception e) {
                throw new RuntimeException("解析CRS出错: " + key, e);
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
     * @param lon WGS84经度（度），用于选择分带
     * @return 逆向投影变换（Gauss→WGS84）
     */
    private MathTransform getTxGaussToWgsByLon(double lon) throws Exception {
        int zone = gkZoneFromLon(lon);
        String key = "G2W:" + zone;
        MathTransform cached = config.txCache.get(key);
        if (cached != null)
            return cached;
        CoordinateReferenceSystem src = getGaussKrugerCRSByLon(lon);
        CoordinateReferenceSystem tgt = getCachedCRS(config.WGS84);
        MathTransform mt = CRS.findMathTransform(src, tgt, true);
        config.txCache.put(key, mt);
        return mt;
    }

    /**
     * 校验 WGS84 坐标几何的经纬度范围是否有效。
     * 经度在 [-180, 180]，纬度在 [-90, 90]；排除 NaN/Inf。
     *
     * @param g WGS84坐标系下的几何
     * @return 坐标范围有效返回 true，否则 false
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
     * 计算两点间速度（km/h）。
     * 使用哈弗辛球面距离与时间差计算瞬时速度；当时间差小于等于0或输入缺失时返回正无穷。
     *
     * @param a 起点（含经纬度与时间）
     * @param b 终点（含经纬度与时间）
     * @return 速度值（km/h）；无法计算时返回 Double.POSITIVE_INFINITY
     */
    private double speedKmh(TrackPoint a, TrackPoint b) {
        if (a == null || b == null || a.getTime() == null || b.getTime() == null) {
            return Double.POSITIVE_INFINITY;
        }
        double meters = haversine(a, b);
        double seconds = (double) java.time.Duration.between(a.getTime(), b.getTime()).toMillis() / 1000.0;
        if (seconds <= 0) {
            return Double.POSITIVE_INFINITY;
        }
        return (meters / seconds) * 3.6; // m/s → km/h
    }

    /**
     * 速度预过滤：保留相邻点瞬时速度落在开区间 (minSpeedKmh, maxSpeedKmh) 的点。
     * 首点始终保留，保持时间升序；边界速度值将被剔除。
     *
     * @param seg         输入轨迹点列表（需包含经纬度与时间，建议已按时间升序）
     * @param minSpeedKmh 最小速度阈值（km/h），严格大于该值才保留
     * @param maxSpeedKmh 最大速度阈值（km/h），严格小于该值才保留
     * @return 过滤后的点列表（时间升序）；当输入少于2个点时返回原列表拷贝
     */
    private List<TrackPoint> filterBySpeedRange(List<TrackPoint> seg, double minSpeedKmh, double maxSpeedKmh) {
        if (seg == null || seg.size() <= 1) {
            return seg == null ? Collections.emptyList() : new ArrayList<>(seg);
        }
        List<TrackPoint> res = new ArrayList<>(seg.size());
        res.add(seg.get(0));
        for (int i = 1; i < seg.size(); i++) {
            TrackPoint prev = seg.get(i - 1);
            TrackPoint curr = seg.get(i);
            double v = speedKmh(prev, curr);
            // 删除小于等于 min 或大于等于 max 的点（删掉 curr）
            if (v > minSpeedKmh && v < maxSpeedKmh) {
                res.add(curr);
            }
        }
        return res;
    }

    /**
     * 线缓冲构建轮廓：在高斯-克吕格米制投影下将轨迹组装为折线，
     * 进行道格拉斯-普克简化后按给定宽度缓冲并合并，最终回转为 WGS84。
     *
     * @param seg    轨迹点列表（WGS84，经纬度与时间有效，建议已按时间升序）
     * @param widthM 总宽度（米），用于缓冲半径计算
     * @return 轮廓几何（WGS84；可能为 MultiPolygon）；输入为空时返回空几何
     * @throws Exception 投影转换、缓冲或几何合并过程中可能抛出的异常
     */
    private Geometry buildOutlineByLineBuffers(List<TrackPoint> seg, double widthM) throws Exception {
        long startTime = System.currentTimeMillis();
        log.trace("[buildOutlineByLineBuffers] 开始处理 {} 个轨迹点，线缓冲宽度: {} 米", seg.size(), widthM);

        if (seg == null || seg.size() < 2) {
            return config.geometryFactory.createGeometryCollection(null);
        }

        double originLon = seg.get(0).getLon();
        MathTransform txWgsToGk = getTxWgsToGaussByLon(originLon);

        // 将点转换到米制坐标并按距离切段
        double breakDist = Math.max(widthM, widthM * config.LINE_BREAK_FACTOR);
        List<List<Coordinate>> lines = new ArrayList<>();
        List<Coordinate> current = new ArrayList<>();
        // 新增：跟踪连续航向变化与累计直线长度
        List<Double> headingWindow = new ArrayList<>();
        double straightRunLen = 0.0;

        long convertTime = 0;
        for (int i = 0; i < seg.size(); i++) {
            TrackPoint p = seg.get(i);
            Coordinate src = new Coordinate(p.getLon(), p.getLat());
            Coordinate c = new Coordinate();
            long t = System.currentTimeMillis();
            JTS.transform(src, c, txWgsToGk);
            convertTime += System.currentTimeMillis() - t;

            if (!current.isEmpty()) {
                Coordinate prev = current.get(current.size() - 1);
                double dx = c.x - prev.x;
                double dy = c.y - prev.y;
                double dist = Math.hypot(dx, dy);
                if (dist > breakDist) {
                    if (current.size() >= 2) {
                        lines.add(current);
                    }
                    current = new ArrayList<>();
                    // 新增：断段后重置直线跑动状态
                    headingWindow.clear();
                    straightRunLen = 0.0;
                }
                // 新增：计算航向并进行“长直线跑动”断段判断
                double hdeg = Math.toDegrees(Math.atan2(dy, dx));
                if (hdeg < 0)
                    hdeg += 360.0;
                headingWindow.add(hdeg);
                if (headingWindow.size() > config.ROAD_STRAIGHT_WINDOW_POINTS) {
                    headingWindow.remove(0);
                }
                straightRunLen += dist;
                if (config.ENABLE_ROAD_STRAIGHT_BREAK && headingWindow.size() >= config.ROAD_STRAIGHT_WINDOW_POINTS) {
                    double maxAbsDiff = 0.0;
                    for (int k = 1; k < headingWindow.size(); k++) {
                        double d = headingWindow.get(k) - headingWindow.get(k - 1);
                        while (d > 180)
                            d -= 360;
                        while (d < -180)
                            d += 360;
                        double ad = Math.abs(d);
                        if (ad > maxAbsDiff)
                            maxAbsDiff = ad;
                    }
                    double straightLenThreshold = Math.max(widthM, widthM * config.ROAD_BREAK_DIST_FACTOR);
                    if (maxAbsDiff < config.ROAD_STRAIGHT_DIFF_THRESHOLD_DEG && straightRunLen > straightLenThreshold) {
                        if (current.size() >= 2) {
                            lines.add(current);
                        }
                        current = new ArrayList<>();
                        headingWindow.clear();
                        straightRunLen = 0.0;
                    }
                }
            }
            current.add(c);
        }
        if (current.size() >= 2) {
            lines.add(current);
        }

        // 构造折线，简化，然后缓冲
        long simplifyTime = 0;
        long bufferTime = 0;
        List<LineString> simplifiedLines = new ArrayList<>();
        BufferParameters params = new BufferParameters();
        params.setQuadrantSegments(config.DEFAULT_BUFFER_QUADRANT);
        params.setEndCapStyle(config.BUFFER_END_CAP_STYLE);
        params.setJoinStyle(config.BUFFER_JOIN_STYLE);
        params.setMitreLimit(config.BUFFER_MITRE_LIMIT);

        double tol = Math.max(0.01, widthM * config.LINE_SIMPLIFY_TOL_FACTOR);

        for (List<Coordinate> coords : lines) {
            if (coords.size() < 2)
                continue;
            Coordinate[] arr = coords.toArray(new Coordinate[0]);
            LineString line = config.geometryFactory.createLineString(arr);
            long ts = System.currentTimeMillis();
            Geometry simplified = DouglasPeuckerSimplifier.simplify(line, tol);
            simplifyTime += System.currentTimeMillis() - ts;

            if (simplified instanceof LineString) {
                LineString ls = (LineString) simplified;
                simplifiedLines.add(ls);
            } else if (simplified instanceof MultiLineString) {
                MultiLineString mls = (MultiLineString) simplified;
                for (int i = 0; i < mls.getNumGeometries(); i++) {
                    LineString ls = (LineString) mls.getGeometryN(i);
                    simplifiedLines.add(ls);
                }
            }
        }

        // 将所有简化线拼为一个 MultiLineString，一次性缓冲
        LineString[] lsArr = simplifiedLines.toArray(new LineString[0]);
        MultiLineString mlsAll = config.geometryFactory.createMultiLineString(lsArr);
        long tb = System.currentTimeMillis();
        Geometry union = BufferOp.bufferOp(mlsAll, widthM, params);
        bufferTime += System.currentTimeMillis() - tb;
        long unionTime = 0;

        // 转回 WGS84
        long backStart = System.currentTimeMillis();
        MathTransform txGkToWgs = getTxGaussToWgsByLon(originLon);
        Geometry result = JTS.transform(union, txGkToWgs);
        long backTime = System.currentTimeMillis() - backStart;

        long endTime = System.currentTimeMillis();
        log.debug("[buildOutlineByLineBuffers] 完成，总计耗时: {}ms (转换:{}ms, 简化:{}ms, 缓冲:{}ms, 合并:{}ms, 回转:{}ms)",
                endTime - startTime, convertTime, simplifyTime, bufferTime, unionTime, backTime);

        return result;
    }

    /**
     * 凹包开运算（保留缝隙/洞）：在米制投影下对轮廓施加轻微开运算，不做闭运算或外缘平滑，
     * 以保留细小缝隙与内部洞结构；处理后回转为 WGS84。
     *
     * @param outlineWgs 输入轮廓（WGS84），允许为 Polygon 或 MultiPolygon
     * @param widthM     总宽度（米），用于推导开运算半径
     * @param originLon  原点经度，用于选取高斯-克吕格投影带（通常取轮廓中心经度）
     * @return 开运算后的几何（WGS84）；输入为空时原样返回
     * @throws Exception 投影转换或几何运算失败时抛出异常
     */
    private Geometry concaveSmoothPreserveHoles(Geometry outlineWgs, double widthM, double originLon) throws Exception {
        if (outlineWgs == null || outlineWgs.isEmpty())
            return outlineWgs;
        long t0 = System.currentTimeMillis();
        MathTransform txWgsToGk = getTxWgsToGaussByLon(originLon);
        MathTransform txGkToWgs = getTxGaussToWgsByLon(originLon);

        double r = Math.max(0.5, widthM * config.CONCAVE_SMOOTH_RADIUS_FACTOR);
        BufferParameters bp = new BufferParameters();
        bp.setQuadrantSegments(config.DEFAULT_BUFFER_QUADRANT);
        bp.setEndCapStyle(config.BUFFER_END_CAP_STYLE);
        bp.setJoinStyle(config.BUFFER_JOIN_STYLE);
        bp.setMitreLimit(config.BUFFER_MITRE_LIMIT);

        Geometry gk = JTS.transform(outlineWgs, txWgsToGk);
        // 改为“开运算”：先负缓冲（收缩以消除窄桥/毛刺），再正缓冲（拉回），更好地保留缝隙与洞
        Geometry opened1 = BufferOp.bufferOp(gk, -r, bp);
        Geometry opened2 = BufferOp.bufferOp(opened1, r, bp);
        Geometry wgs = JTS.transform(opened2, txGkToWgs);
        long t1 = System.currentTimeMillis();
        log.debug("[concaveSmooth:opening] 半径={}m 输入类型={} 输出类型={} 耗时={}ms", r, outlineWgs.getGeometryType(),
                wgs.getGeometryType(), (t1 - t0));
        return wgs;
    }

    /**
     * 检查几何是否包含有效坐标。
     * 要求坐标数至少为1，并且不含 NaN/Inf。
     *
     * @param g 几何
     * @return 坐标有效返回 true，否则 false
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
     * 保留面积最大的前 N 个区块
     * 适用于 `Polygon` 或 `MultiPolygon`：按 JTS 计算的平方米面积排序并截取前 `maxCount` 个。
     *
     * @param geometry 输入几何（`Polygon`/`MultiPolygon`），为空或非面类型时直接返回原值
     * @param maxCount 保留的最大区块数量（<=0 时返回原值）
     * @return 面积前 N 的几何；当 N 小于区块数时可能返回 `MultiPolygon`
     */
    private Geometry keepLargestPolygons(Geometry geometry, int maxCount) {
        if (geometry == null)
            return null;
        if (maxCount <= 0)
            return geometry;

        GeometryFactory gf = geometry.getFactory();
        List<Polygon> polys = new ArrayList<>();

        // 展开并收集所有 Polygon
        for (int i = 0; i < geometry.getNumGeometries(); i++) {
            Geometry g = geometry.getGeometryN(i);
            if (g instanceof Polygon) {
                polys.add((Polygon) g);
            } else if (g instanceof MultiPolygon) {
                MultiPolygon mp = (MultiPolygon) g;
                for (int j = 0; j < mp.getNumGeometries(); j++) {
                    polys.add((Polygon) mp.getGeometryN(j));
                }
            }
            // 其他类型忽略
        }
        if (polys.isEmpty())
            return geometry;

        // 按面积倒序
        polys.sort(Comparator.comparingDouble(Polygon::getArea).reversed());

        int limit = Math.min(maxCount, polys.size());
        if (limit == 1) {
            return polys.get(0);
        }
        Polygon[] top = polys.subList(0, limit).toArray(new Polygon[0]);
        return gf.createMultiPolygon(top);
    }

    /**
     * 按面积阈值过滤 MultiPolygon 中的小区块（亩）。
     * 单个 Polygon 不进行过滤；MultiPolygon 中逐个子面计算球面面积（与 Turf.js 对齐）并换算为亩，移除小于 `minMu`
     * 的区块。
     *
     * @param geometry 输入几何（`Polygon`/`MultiPolygon`）
     * @param minMu    最小亩数阈值（>0），小于等于 0 时返回原值
     * @return 过滤后的几何；若全部被移除则返回空几何
     */
    private Geometry removeSmallMuPolygons(Geometry geometry, double minMu) {
        if (geometry == null)
            return null;
        GeometryFactory gf = geometry.getFactory();
        if (geometry instanceof Polygon) {
            Polygon poly = (Polygon) geometry;
            // 单个 Polygon 不进行最小亩数过滤，直接保留
            return poly;
        } else if (geometry instanceof MultiPolygon) {
            MultiPolygon mp = (MultiPolygon) geometry;
            List<Polygon> kept = new ArrayList<>();
            for (int i = 0; i < mp.getNumGeometries(); i++) {
                Polygon p = (Polygon) mp.getGeometryN(i);
                double mu = calcMu(p);
                if (mu >= minMu) {
                    kept.add(p);
                }
            }
            if (kept.isEmpty()) {
                return gf.createMultiPolygon(new Polygon[0]);
            }
            if (kept.size() == 1) {
                return kept.get(0);
            }
            return gf.createMultiPolygon(kept.toArray(new Polygon[0]));
        }
        return geometry;
    }

    /**
     * 外缘细长裁剪：仅基于外环进行轻微开运算以移除边缘细长条（在米制投影下执行）。
     * 仅在 `splitRoad` 中可能启用；`getOutline` 不使用该裁剪。
     *
     * @param g       输入几何（WGS84）
     * @param radiusM 开运算半径（米）
     * @return 裁剪后的几何（WGS84）
     */
    private Geometry trimOuterThinStrips(Geometry g, double radiusM) {
        if (g == null || g.isEmpty())
            return g;
        try {
            // 根据几何包络的中心经度选择高斯-克吕格分带
            Envelope env = g.getEnvelopeInternal();
            double cenLon = env.getMinX() + env.getWidth() / 2.0;
            MathTransform txWgsToGk = getTxWgsToGaussByLon(cenLon);
            MathTransform txGkToWgs = getTxGaussToWgsByLon(cenLon);

            Geometry gk = JTS.transform(g, txWgsToGk);
            GeometryFactory gf = gk.getFactory();

            // 收集所有外环组成统一外壳几何，避免逐面两次缓冲的高成本
            List<Polygon> shells = new ArrayList<>();
            if (gk instanceof Polygon) {
                Polygon p = (Polygon) gk;
                shells.add(gf.createPolygon(p.getExteriorRing().getCoordinates()));
            } else if (gk instanceof MultiPolygon) {
                MultiPolygon mp = (MultiPolygon) gk;
                for (int i = 0; i < mp.getNumGeometries(); i++) {
                    Polygon p = (Polygon) mp.getGeometryN(i);
                    shells.add(gf.createPolygon(p.getExteriorRing().getCoordinates()));
                }
            } else {
                return g;
            }
            if (shells.isEmpty())
                return g;
            Geometry shellsGeomGk = shells.size() == 1
                    ? shells.get(0)
                    : gf.createMultiPolygon(shells.toArray(new Polygon[0]));

            // 可选：对外环集合做拓扑保留简化以减小缓冲成本（容差不超过1.0m，随半径调整）
            double tolThin = Math.min(1.0, radiusM * 0.5);
            Geometry shellsSimplified = TopologyPreservingSimplifier.simplify(shellsGeomGk, tolThin);
            if (shellsSimplified == null || shellsSimplified.isEmpty())
                shellsSimplified = shellsGeomGk;

            // 使用配置化缓冲参数，统一执行一次开运算（负缓冲再正缓冲）
            BufferParameters params = new BufferParameters(
                    config.DEFAULT_BUFFER_QUADRANT,
                    config.BUFFER_END_CAP_STYLE,
                    config.BUFFER_JOIN_STYLE,
                    config.BUFFER_MITRE_LIMIT);
            Geometry eroded = BufferOp.bufferOp(shellsSimplified, -radiusM, params);
            if (eroded == null || eroded.isEmpty()) {
                // 半径过大导致核心消失，保持原面
                return g;
            }
            Geometry reopened = BufferOp.bufferOp(eroded, radiusM, params);

            // 使用 PreparedGeometry，逐面裁剪以降低一次性大规模交集的开销
            PreparedGeometry prepReopened = PreparedGeometryFactory.prepare(reopened);
            java.util.List<Geometry> trimmedParts = new java.util.ArrayList<>();
            if (gk instanceof Polygon) {
                Polygon p = (Polygon) gk;
                if (prepReopened.intersects(p)) {
                    Geometry t = p.intersection(reopened);
                    if (t != null && !t.isEmpty())
                        trimmedParts.add(t);
                }
            } else if (gk instanceof MultiPolygon) {
                MultiPolygon mp2 = (MultiPolygon) gk;
                for (int i = 0; i < mp2.getNumGeometries(); i++) {
                    Polygon p = (Polygon) mp2.getGeometryN(i);
                    if (prepReopened.intersects(p)) {
                        Geometry t = p.intersection(reopened);
                        if (t != null && !t.isEmpty())
                            trimmedParts.add(t);
                    }
                }
            }
            Geometry trimmedGk = trimmedParts.isEmpty() ? gk
                    : UnaryUnionOp.union(trimmedParts);
            return JTS.transform(trimmedGk, txGkToWgs);
        } catch (Exception ex) {
            log.warn("[trimOuterThin] 处理失败: {}", ex.getMessage());
            return g;
        }
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
                return "GEOMETRYCOLLECTION EMPTY";
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
            return g != null ? g.toText() : "GEOMETRYCOLLECTION EMPTY";
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
        Geometry inter = g1.intersection(g2);
        if (inter == null || inter.isEmpty()) {
            return new WktIntersectionResult(null, 0.0);
        }
        String wkt = toWkt(inter);
        double mu = calcMu(inter);
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
     * @throws Exception 坐标转换或几何运算异常
     */
    public OutlinePart getOutline(List<TrackPoint> seg, double totalWidthM) throws Exception {
        long t0 = System.currentTimeMillis();
        if (seg == null) {
            throw new IllegalArgumentException("轨迹点列表不能为空");
        }
        if (totalWidthM < 0) {
            throw new IllegalArgumentException("总宽度必须为非负数");
        }
        double widthM = totalWidthM / 2.0;

        // 仅过滤异常点（越界与(0,0)），按时间升序；不做切割与删除
        List<TrackPoint> points = filterAndSortTrackPoints(seg);
        log.trace("[getOutline] 轨迹点过滤完成，原始点数: {}, 过滤后点数: {}", seg.size(), points.size());
        if (points.size() < 3) {
            throw new IllegalArgumentException("轨迹段至少需要3个有效点");
        }

        // 计算起止时间（忽略空时间）
        LocalDateTime startTime = null;
        LocalDateTime endTime = null;
        for (TrackPoint p : points) {
            if (p.getTime() != null) {
                if (startTime == null)
                    startTime = p.getTime();
                endTime = p.getTime();
            }
        }

        // 使用线缓冲轮廓构建（WGS84坐标系），保持过滤异常点与不拆分道路，只返回单个Polygon
        long tBuildStart = System.currentTimeMillis();
        Geometry outline = buildOutlineByLineBuffers(points, widthM);
        long tBuildEnd = System.currentTimeMillis();
        log.trace("[getOutline] 轮廓构建完成（线缓冲），耗时: {}ms", (tBuildEnd - tBuildStart));

        // 基于首点经度构建高斯-克吕格米制投影转换（形态学处理）
        double originLonGk = points.get(0).getLon();
        MathTransform txWgsToGkOutline = getTxWgsToGaussByLon(originLonGk);
        MathTransform txGkToWgsOutline = getTxGaussToWgsByLon(originLonGk);

        // 若为 MultiPolygon，则转换为单个 Polygon（凹包）：基于 alpha-like 形态操作
        if (outline instanceof MultiPolygon) {
            try {
                // 在米制下进行形态学处理以获得凹包
                Geometry unioned = UnaryUnionOp.union(outline);
                Geometry proj = JTS.transform(unioned, txWgsToGkOutline);
                BufferParameters bp = new BufferParameters();
                bp.setQuadrantSegments(config.DEFAULT_BUFFER_QUADRANT);
                bp.setEndCapStyle(BufferParameters.CAP_FLAT);
                bp.setJoinStyle(BufferParameters.JOIN_BEVEL);
                bp.setMitreLimit(2.0);

                // 不进行凸包/闭运算，不裁剪不平滑；仅做轻微“开运算”以保留缝隙
                double openR = Math.max(0.2, widthM * config.GAP_ENHANCE_ERODE_FACTOR);
                Geometry opened = BufferOp.bufferOp(proj, -openR, bp);
                opened = BufferOp.bufferOp(opened, openR, bp);
                Geometry openedWgs = JTS.transform(opened, txGkToWgsOutline);

                if (openedWgs instanceof Polygon) {
                    outline = (Polygon) openedWgs;
                } else if (openedWgs instanceof MultiPolygon) {
                    // 保留缝隙与洞，不做闭运算；选择面积最大的外轮廓
                    MultiPolygon mp = (MultiPolygon) openedWgs;
                    int maxIdx = 0;
                    double maxArea = -1;
                    for (int i = 0; i < mp.getNumGeometries(); i++) {
                        Polygon p = (Polygon) mp.getGeometryN(i);
                        double a = p.getArea();
                        if (a > maxArea) {
                            maxArea = a;
                            maxIdx = i;
                        }
                    }
                    outline = (Polygon) mp.getGeometryN(maxIdx);
                }
                log.trace("[getOutline] 已将MultiPolygon转换为单个Polygon（凹包，保留缝隙，未裁剪/未平滑）");
            } catch (Exception ex) {
                log.warn("[getOutline] 凹包处理失败，保留最大面积外轮廓: {}", ex.getMessage());
                if (outline instanceof MultiPolygon) {
                    MultiPolygon mp = (MultiPolygon) outline;
                    int maxIdx = 0;
                    double maxArea = -1;
                    for (int i = 0; i < mp.getNumGeometries(); i++) {
                        Polygon p = (Polygon) mp.getGeometryN(i);
                        double a = p.getArea();
                        if (a > maxArea) {
                            maxArea = a;
                            maxIdx = i;
                        }
                    }
                    outline = (Polygon) mp.getGeometryN(maxIdx);
                } else if (outline instanceof Polygon) {
                    // 已是单面，直接保留
                } else {
                    Geometry env = outline.getEnvelope();
                    outline = (env instanceof Polygon)
                            ? (Polygon) env
                            : (Polygon) ((MultiPolygon) env).getGeometryN(0);
                }
            }
        }

        // 计算亩数与WKT（不删除任何区块）
        double mu = calcMu(outline);
        String wkt = toWkt(outline);

        // 记录结果类型
        log.debug("[getOutline] 返回结果 type=Polygon");

        long t1 = System.currentTimeMillis();
        log.debug("[getOutline] 总耗时={}ms", (t1 - t0));
        return new OutlinePart(outline, startTime, endTime, mu, wkt, points).setTotalWidthM(totalWidthM);
    }

    /**
     * 基于轨迹点与总宽度生成轮廓，并按面积保留 Top-N 最大区块。
     *
     * @param seg         轨迹点列表（WGS84），至少 3 个有效点
     * @param totalWidthM 道路总宽度（米），左右合计，必须非负
     * @return 轮廓几何，可能为 `Polygon` 或 `MultiPolygon`（裁剪后最多保留前 N 个）
     * @throws Exception 坐标转换或几何运算异常
     */
    public SplitRoadResult splitRoad(List<TrackPoint> seg, double totalWidthM) throws Exception {
        return splitRoad(seg, totalWidthM, config.DEFAULT_MAX_OUTLINE_SEGMENTS);
    }

    /**
     * 基于轨迹点与总宽度生成轮廓，并按面积保留不超过 maxSegments 的最大区块。
     * 
     * @param seg         轨迹点列表（WGS84），至少 3 个有效点
     * @param totalWidthM 作业机具总宽度（米），左右合计，必须非负
     * @param maxSegments 返回的区块数量上限；为 null 或 <=0 时取默认值 N=10
     * @return 轮廓几何，可能为 `Polygon` 或 `MultiPolygon`（最多保留 `maxSegments` 个）
     * @throws Exception 坐标转换或几何运算异常
     */
    public SplitRoadResult splitRoad(List<TrackPoint> seg, double totalWidthM, Integer maxSegments) throws Exception {
        long t0 = System.currentTimeMillis();
        double widthM = totalWidthM / 2.0;
        log.debug("[splitRoad] 开始 参数 原始点数={} 总宽度={}m 返回上限={}", seg.size(), totalWidthM,
                (maxSegments == null || maxSegments <= 0) ? config.DEFAULT_MAX_OUTLINE_SEGMENTS
                        : maxSegments.intValue());

        // 1) 过滤异常点（越界与(0,0)），按时间升序
        log.debug("[splitRoad] 异常点过滤 开始 参数 输入点数={}", seg.size());
        long tFilterStart = System.currentTimeMillis();
        List<TrackPoint> sortedSeg = filterAndSortTrackPoints(seg);
        long tFilterEnd = System.currentTimeMillis();
        log.debug("[splitRoad] 异常点过滤 结束 返回 返回点数={} 耗时={}ms", sortedSeg.size(), (tFilterEnd - tFilterStart));

        // 2) 速度过滤：删除小于等于1km/h或大于等于20km/h的点（使用配置阈值）
        log.debug("[splitRoad] 速度过滤 开始 范围=[{}km/h, {}km/h] 输入点数={}",
                config.MIN_WORK_SPEED_KMH, config.WORK_MAX_SPEED_KMH, sortedSeg.size());
        long tSpeedStart = System.currentTimeMillis();
        List<TrackPoint> workSeg = filterBySpeedRange(sortedSeg, config.MIN_WORK_SPEED_KMH, config.WORK_MAX_SPEED_KMH);
        long tSpeedEnd = System.currentTimeMillis();
        int removed = sortedSeg.size() - workSeg.size();
        log.debug("[splitRoad] 速度过滤 结束 返回 移除点数={} 剩余点数={} 耗时={}ms",
                removed, workSeg.size(), (tSpeedEnd - tSpeedStart));

        // 3) 直接进行线缓冲构建轮廓（不使用点缓冲）
        log.debug("[splitRoad] 构建线缓冲 开始 宽度={}m 输入点数={}", widthM, workSeg.size());
        long tBuildStart = System.currentTimeMillis();
        Geometry outline = buildOutlineByLineBuffers(workSeg, widthM);
        long tBuildEnd = System.currentTimeMillis();
        log.debug("[splitRoad] 构建线缓冲 结束 返回 type={} parts={} 耗时={}ms", outline.getGeometryType(),
                (outline instanceof MultiPolygon) ? outline.getNumGeometries() : 1, (tBuildEnd - tBuildStart));
        int partsOutline = (outline instanceof MultiPolygon) ? outline.getNumGeometries() : 1;

        // 3.1) 凹包平滑（保留镂空/洞），轻微开运算以平滑外缘并扩大缝隙
        if (config.ENABLE_CONCAVE_SMOOTHING && workSeg != null && !workSeg.isEmpty()) {
            long tSmoothStart = System.currentTimeMillis();
            log.debug("[splitRoad] 凹包平滑 开始 半径系数={} 宽度={}m 输入类型={}", config.CONCAVE_SMOOTH_RADIUS_FACTOR, widthM,
                    outline.getGeometryType());
            double originLon = workSeg.get(0).getLon();
            outline = concaveSmoothPreserveHoles(outline, widthM, originLon);
            long tSmoothEnd = System.currentTimeMillis();
            log.debug("[splitRoad] 凹包平滑 结束 返回 输出类型={} 耗时={}ms", outline.getGeometryType(),
                    (tSmoothEnd - tSmoothStart));
        }

        // 3.2) 外缘细长条裁剪（避免道路窄条并入田块外缘）
        if (config.ENABLE_OUTER_THIN_TRIM) {
            long tThinStart = System.currentTimeMillis();
            double radiusM = Math.max(0.5, widthM * config.THIN_TRIM_RADIUS_FACTOR);
            int partsBeforeThin = (outline instanceof MultiPolygon) ? outline.getNumGeometries() : 1;
            log.debug("[splitRoad] 外缘细长裁剪 开始 半径={}m 输入类型={} 输入部分数={}", radiusM, outline.getGeometryType(),
                    partsBeforeThin);
            Geometry outlineAfterThin = trimOuterThinStrips(outline, radiusM);
            int partsAfterThin = (outlineAfterThin instanceof MultiPolygon) ? outlineAfterThin.getNumGeometries() : 1;
            int removedThin = Math.max(0, partsBeforeThin - partsAfterThin);
            outline = outlineAfterThin;
            long tThinEnd = System.currentTimeMillis();
            log.debug("[splitRoad] 外缘细长裁剪 结束 返回 输出类型={} 移除部分数={} 剩余部分数={} 耗时={}ms", outline.getGeometryType(),
                    removedThin, partsAfterThin,
                    (tThinEnd - tThinStart));
        }

        // 3.3) 缝隙扩大轻微蚀刻（小幅整体收缩以凸显镂空与缝隙）
        if (config.ENABLE_GAP_ENHANCE_ERODE && workSeg != null && !workSeg.isEmpty()) {
            long tErodeStart = System.currentTimeMillis();
            double originLon = workSeg.get(0).getLon();
            double erodeR = Math.max(0.2, widthM * config.GAP_ENHANCE_ERODE_FACTOR);
            log.debug("[splitRoad] 缝隙扩大蚀刻 开始 半径={}m 因子={} 输入类型={}", erodeR, config.GAP_ENHANCE_ERODE_FACTOR,
                    outline.getGeometryType());
            MathTransform txWgsToGk = getTxWgsToGaussByLon(originLon);
            MathTransform txGkToWgs = getTxGaussToWgsByLon(originLon);
            Geometry outlineGk = JTS.transform(outline, txWgsToGk);
            BufferParameters eParams = new BufferParameters();
            eParams.setQuadrantSegments(config.DEFAULT_BUFFER_QUADRANT);
            eParams.setEndCapStyle(config.BUFFER_END_CAP_STYLE);
            eParams.setJoinStyle(config.BUFFER_JOIN_STYLE);
            eParams.setMitreLimit(config.BUFFER_MITRE_LIMIT);
            Geometry erodedGk = BufferOp.bufferOp(outlineGk, -erodeR, eParams);
            outline = JTS.transform(erodedGk, txGkToWgs);
            long tErodeEnd = System.currentTimeMillis();
            log.debug("[splitRoad] 缝隙扩大蚀刻 结束 返回 输出类型={} 耗时={}ms", outline.getGeometryType(),
                    (tErodeEnd - tErodeStart));
        }

        // 4) 保留最大区块并按最小亩数过滤
        int limit = (maxSegments == null || maxSegments <= 0) ? config.DEFAULT_MAX_OUTLINE_SEGMENTS
                : maxSegments.intValue();
        log.debug("[splitRoad] 输入参数 原始点数={} 过滤后点数={} 总宽度={}m 返回上限={}", seg.size(), workSeg.size(), totalWidthM, limit);
        log.debug("[splitRoad] 保留最大区块 开始 上限={} 输入类型={} parts={}", limit, outline.getGeometryType(), partsOutline);
        long tKeepStart = System.currentTimeMillis();
        Geometry kept = keepLargestPolygons(outline, limit);
        long tKeepEnd = System.currentTimeMillis();
        int partsAfterKeep = (kept instanceof MultiPolygon) ? kept.getNumGeometries() : 1;
        log.debug("[splitRoad] 保留最大区块 结束 返回 输出部分数={} 耗时={}ms", partsAfterKeep, (tKeepEnd - tKeepStart));

        log.trace("[splitRoad] 动态最小亩数阈值 minMuDynamic={} (partsAfterKeep={})", config.MIN_MU_DYNAMIC_THRESHOLD_MU,
                partsAfterKeep);
        log.debug("[splitRoad] 小面积过滤 开始 阈值(亩)={} 输入部分数={}", config.MIN_MU_DYNAMIC_THRESHOLD_MU, partsAfterKeep);
        long tAreaStart = System.currentTimeMillis();
        Geometry trimmed = removeSmallMuPolygons(kept, config.MIN_MU_DYNAMIC_THRESHOLD_MU);
        long tAreaEnd = System.currentTimeMillis();
        int partsAfterAreaFilter = (trimmed instanceof MultiPolygon) ? trimmed.getNumGeometries() : 1;
        int removedByArea = partsAfterKeep - partsAfterAreaFilter;
        log.debug("[splitRoad] 小面积过滤 结束 返回 移除区块数={} 剩余部分数={} 耗时={}ms", removedByArea, partsAfterAreaFilter,
                (tAreaEnd - tAreaStart));
        log.trace("[splitRoad] 裁剪完成 type={} parts={} 移除面积={} 耗时={}ms",
                trimmed.getGeometryType(), (trimmed instanceof MultiPolygon) ? trimmed.getNumGeometries() : 1,
                removedByArea, ((tKeepEnd - tKeepStart) + (tAreaEnd - tAreaStart)));

        // 5) 返回日志说明
        if (((trimmed instanceof MultiPolygon) ? trimmed.getNumGeometries() : 1) < limit) {
            String reason;
            if (((trimmed instanceof MultiPolygon) ? trimmed.getNumGeometries() : 1) != partsAfterKeep) {
                reason = "小面积过滤";
            } else if (partsOutline < limit) {
                reason = "原始区块数小于上限";
            } else {
                reason = "面积裁剪后区块数不超过上限";
            }
            log.debug("[splitRoad] 返回结果 结束 返回 type={} parts={} 原因={} 小面积移除={}", trimmed.getGeometryType(),
                    (trimmed instanceof MultiPolygon) ? trimmed.getNumGeometries() : 1, reason, removedByArea);
        } else {
            log.debug("[splitRoad] 返回结果 结束 返回 type={} parts={}", trimmed.getGeometryType(),
                    (trimmed instanceof MultiPolygon) ? trimmed.getNumGeometries() : 1);
        }

        // 6) 统计轮廓内点的会话段（与轮廓构建一致的过滤点集）
        List<OutlinePart> parts = new ArrayList<>();
        if (trimmed instanceof Polygon) {
            log.trace("[splitRoad] 进入Polygon分支");
            Polygon poly = (Polygon) trimmed;
            PreparedGeometry preparedPoly = PreparedGeometryFactory.prepare(poly);
            Envelope env = poly.getEnvelopeInternal();
            long tWithinStart = System.currentTimeMillis();
            List<List<TrackPoint>> sessions = new ArrayList<>();
            List<TrackPoint> current = null;
            for (TrackPoint p : workSeg) {
                Coordinate c = new Coordinate(p.getLon(), p.getLat());
                if (!env.contains(c)) {
                    if (current != null) {
                        sessions.add(current);
                        current = null;
                    }
                    continue;
                }
                Geometry point = config.geometryFactory.createPoint(c);
                boolean inside = preparedPoly.contains(point);
                if (inside) {
                    if (current == null)
                        current = new ArrayList<>();
                    current.add(p);
                } else {
                    if (current != null) {
                        sessions.add(current);
                        current = null;
                    }
                }
            }
            if (current != null) {
                sessions.add(current);
            }
            long tWithinEnd = System.currentTimeMillis();
            log.trace("[splitRoad] Polygon 会话段统计完成 sessions={} 耗时={}ms", sessions.size(), (tWithinEnd - tWithinStart));

            List<TrackPoint> inPoly = sessions.stream().flatMap(List::stream).collect(Collectors.toList());
            LocalDateTime start = inPoly.isEmpty() ? null : inPoly.get(0).getTime();
            LocalDateTime end = inPoly.isEmpty() ? null : inPoly.get(inPoly.size() - 1).getTime();
            log.trace("[splitRoad] Polygon 跨度(按轮廓内点) 时间范围 start={} end={} inPoly={} sessions={}", start, end,
                    inPoly.size(), sessions.size());

            long tMuStart = System.currentTimeMillis();
            double mu = calcMu(poly);
            long tMuEnd = System.currentTimeMillis();
            log.trace("[splitRoad] Polygon mu计算完成 值={} 耗时={}ms", mu, (tMuEnd - tMuStart));

            long tWktStart = System.currentTimeMillis();
            String wkt = toWkt(poly);
            long tWktEnd = System.currentTimeMillis();
            log.trace("[splitRoad] Polygon WKT生成完成 长度={} 耗时={}ms", wkt.length(), (tWktEnd - tWktStart));

            parts.add(new OutlinePart(poly, start, end, mu, wkt, inPoly).setTotalWidthM(totalWidthM));
            log.trace("[splitRoad] parts构建完成 size={}", parts.size());
        } else if (trimmed instanceof MultiPolygon) {
            MultiPolygon mp = (MultiPolygon) trimmed;
            log.trace("[splitRoad] 进入MultiPolygon分支 分区数={}", mp.getNumGeometries());
            long partsBuildStart = System.currentTimeMillis();
            for (int i = 0; i < mp.getNumGeometries(); i++) {
                Polygon poly = (Polygon) mp.getGeometryN(i);
                PreparedGeometry preparedPoly = PreparedGeometryFactory.prepare(poly);
                Envelope env = poly.getEnvelopeInternal();
                long tWithinStart = System.currentTimeMillis();
                List<List<TrackPoint>> sessions = new ArrayList<>();
                List<TrackPoint> current = null;
                for (TrackPoint p : workSeg) {
                    Coordinate c = new Coordinate(p.getLon(), p.getLat());
                    if (!env.contains(c)) {
                        if (current != null) {
                            sessions.add(current);
                            current = null;
                        }
                        continue;
                    }
                    Geometry point = config.geometryFactory.createPoint(c);
                    boolean inside = preparedPoly.contains(point);
                    if (inside) {
                        if (current == null)
                            current = new ArrayList<>();
                        current.add(p);
                    } else {
                        if (current != null) {
                            sessions.add(current);
                            current = null;
                        }
                    }
                }
                if (current != null) {
                    sessions.add(current);
                }
                long tWithinEnd = System.currentTimeMillis();
                log.trace("[splitRoad] Part#{} 会话段统计完成 sessions={} 耗时={}ms", i, sessions.size(),
                        (tWithinEnd - tWithinStart));

                List<TrackPoint> inPoly = sessions.stream().flatMap(List::stream).collect(Collectors.toList());
                LocalDateTime start = inPoly.isEmpty() ? null : inPoly.get(0).getTime();
                LocalDateTime end = inPoly.isEmpty() ? null : inPoly.get(inPoly.size() - 1).getTime();
                log.trace("[splitRoad] Part#{} 跨度(按轮廓内点) 时间范围 start={} end={} (sessions={})", i, start, end,
                        sessions.size());

                long tMuStart = System.currentTimeMillis();
                double mu = calcMu(poly);
                long tMuEnd = System.currentTimeMillis();
                log.trace("[splitRoad] Part#{} mu计算完成 值={} 耗时={}ms", i, mu, (tMuEnd - tMuStart));

                long tWktStart = System.currentTimeMillis();
                String wkt = toWkt(poly);
                long tWktEnd = System.currentTimeMillis();
                log.trace("[splitRoad] Part#{} WKT生成完成 长度={} 耗时={}ms", i, wkt.length(), (tWktEnd - tWktStart));

                parts.add(new OutlinePart(poly, start, end, mu, wkt, inPoly).setTotalWidthM(totalWidthM));
                if ((i + 1) % 10 == 0 || i == mp.getNumGeometries() - 1) {
                    log.trace("[splitRoad] MultiPolygon 进度 {}/{}", (i + 1), mp.getNumGeometries());
                }
            }
            long partsBuildEnd = System.currentTimeMillis();
            log.trace("[splitRoad] MultiPolygon parts构建完成 size={} 总耗时={}ms", parts.size(),
                    (partsBuildEnd - partsBuildStart));
        }

        long tOutlineWktStart = System.currentTimeMillis();
        String outlineWkt = toWkt(trimmed);
        long tOutlineWktEnd = System.currentTimeMillis();
        log.trace("[splitRoad] Outline WKT生成完成 长度={} 耗时={}ms", outlineWkt.length(),
                (tOutlineWktEnd - tOutlineWktStart));

        long tTotalEnd = System.currentTimeMillis();
        log.debug("[splitRoad] 结束 返回 parts={} 总耗时={}ms",
                (trimmed instanceof MultiPolygon) ? trimmed.getNumGeometries() : 1, (tTotalEnd - t0));
        return new SplitRoadResult(trimmed, parts, outlineWkt).setTotalWidthM(totalWidthM);
    }

}