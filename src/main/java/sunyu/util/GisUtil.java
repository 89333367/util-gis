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
import org.locationtech.jts.geom.MultiPolygon;
import org.locationtech.jts.geom.Polygon;
import org.locationtech.jts.geom.prep.PreparedGeometry;
import org.locationtech.jts.geom.prep.PreparedGeometryFactory;
import org.locationtech.jts.io.ParseException;
import org.locationtech.jts.io.WKTReader;
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
 * <p>
 * 公共方法一览：
 * </p>
 * <ul>
 * <li>builder()：创建构建器，用于配置并构建 {@code GisUtil}</li>
 * <li>close()：释放资源（实现 {@code AutoCloseable}）</li>
 * <li>toWkt(Geometry)：将几何转换为 WKT（统一 WGS84），含坐标识别与修复</li>
 * <li>fromWkt(String)：解析 WKT（WGS84）为 Geometry 并转换到高斯-克吕格米制坐标</li>
 * <li>haversine(CoordinatePoint, CoordinatePoint)：计算两点大圆距离（米，WGS84）</li>
 * <li>calcMu(Geometry)：计算几何面积的亩数（球面公式）</li>
 * <li>calcMu(String)：解析 WKT 后计算亩数</li>
 * <li>intersection(String, String)：计算两 WKT 相交的几何与面积，返回
 * {@code WktIntersectionResult}</li>
 * <li>intersects(String, String)：判断两 WKT 是否相交（WGS84，支持
 * {@code POLYGON}/{@code MULTIPOLYGON}）</li>
 * <li>equalsWkt(String, String)：判断两 WKT 是否拓扑相等</li>
 * <li>disjoint(String, String)：判断是否脱节</li>
 * <li>touches(String, String)：判断是否接触（边界接触）</li>
 * <li>crosses(String, String)：判断是否交叉</li>
 * <li>within(String, String)：判断 A 是否在 B 内</li>
 * <li>contains(String, String)：判断 A 是否包含 B</li>
 * <li>overlaps(String, String)：判断是否重叠</li>
 * <li>pointInPolygon(CoordinatePoint, String)：判断点是否在多边形内（含边界）</li>
 * <li>splitRoad(List&lt;TrackPoint&gt;, double)：按总宽度对轨迹进行分段并返回结果</li>
 * <li>splitRoad(List&lt;TrackPoint&gt;, double, Integer)：指定最大段数的道路分段</li>
 * </ul>
 * 
 * <p>
 * 概览：
 * </p>
 * <ul>
 * <li>坐标系：默认 WGS84；轮廓构建与形态学在高斯-克吕格米制投影（按首点经度分带）进行，结果回转 WGS84。</li>
 * <li>常量：距离计算使用平均地球半径 R=6371000；面积计算在 ringArea 使用 WGS84 半径 6378137（与 Turf.js
 * 对齐）。</li>
 * <li>轮廓：线缓冲构建（折线简化+缓冲，拐角/端头由 BufferParameters 控制）；splitRoad
 * 可选外缘细长裁剪与缝隙增强蚀刻。</li>
 * <li>边界：{@code pointInPolygon} 边界视为内；几何谓词遵循 JTS 语义。</li>
 * </ul>
 * 
 * <p>
 * 设计要点：
 * </p>
 * <ul>
 * <li>坐标系：内部优先在 WGS84 下处理；涉及形态学与面积计算时使用高斯-克吕格米制投影（6 度分带）。</li>
 * <li>变换缓存：按分带缓存 CRS 与 MathTransform，避免重复构建。</li>
 * <li>轮廓构建：线简化+线缓冲并合并，控制几何规模与性能。</li>
 * <li>清理过滤：按面积（亩数）等规则过滤过小碎片（仅 splitRoad 中使用）。</li>
 * <li>输出：支持 WKT 输出并在必要时进行坐标系识别与修复。</li>
 * </ul>
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

        // 默认轮廓返回的最多多边形数量（TopN）
        private final int DEFAULT_MAX_OUTLINE_SEGMENTS = 10;

        // 圆近似细分：4，提升性能，误差满足道路宽度场景
        private final int DEFAULT_BUFFER_QUADRANT = 2;

        // 作业最大速度阈值（km/h），用于前置速度过滤；可通过 Builder 配置
        private double WORK_MAX_SPEED_KMH = 20.0;
        // 作业最小速度阈值（km/h），用于前置速度过滤；可通过 Builder 配置
        private double MIN_WORK_SPEED_KMH = 1.0;

        // 距离断开阈值因子（倍数*单侧宽度），用于buildOutlineByLineBuffersWithDistanceBreak方法
        private double DISTANCE_BREAK_FACTOR = 1.5;

        // 最小亩数阈值（亩），小于此值的多边形将被过滤掉
        private double MIN_AREA_MU = 0.5;

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
         * 设置距离断开阈值因子（倍数*单侧宽度），用于buildOutlineByLineBuffersWithDistanceBreak方法。
         * 当相邻轨迹点距离超过（单侧宽度*此因子）时，会断开窗口形成独立多边形。
         * 
         * @param factor 距离断开阈值因子，建议值：2.0（默认），可根据需要调整为3.0或其他值
         * @return Builder实例，支持链式调用
         */
        public Builder withDistanceBreakFactor(double factor) {
            if (factor <= 0) {
                throw new IllegalArgumentException("距离断开阈值因子必须大于0");
            }
            config.DISTANCE_BREAK_FACTOR = factor;
            return this;
        }

        /**
         * 设置最小亩数阈值（亩），小于此值的多边形将被过滤掉。
         * 
         * @param minAreaMu 最小亩数阈值，建议值：0.5（默认），可根据需要调整
         * @return Builder实例，支持链式调用
         */
        public Builder withMinAreaMu(double minAreaMu) {
            if (minAreaMu < 0) {
                throw new IllegalArgumentException("最小亩数阈值必须大于等于0");
            }
            config.MIN_AREA_MU = minAreaMu;
            return this;
        }

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
     * 构建线缓冲轮廓，当相邻点距离超过阈值时断开窗口形成独立多边形
     * 
     * @param points       轨迹点列表
     * @param widthM       缓冲宽度（米）
     * @param maxDistanceM 最大允许距离（米），超过此距离断开窗口
     * @return 几何轮廓（Polygon或MultiPolygon）
     * @throws Exception 坐标转换或几何运算异常
     */
    private Geometry buildOutlineByLineBuffersWithDistanceBreak(List<TrackPoint> points, double widthM,
            double maxDistanceM) throws Exception {
        if (points == null || points.size() < 2) {
            return config.geometryFactory.createPolygon();
        }

        // 坐标转换到米制
        List<TrackPoint> meterPoints = new ArrayList<>();
        for (TrackPoint tp : points) {
            double originLon = tp.getLon();
            MathTransform txWgsToGk = getTxWgsToGaussByLon(originLon);
            Coordinate src = new Coordinate(tp.getLon(), tp.getLat());
            Coordinate c = new Coordinate();
            JTS.transform(src, c, txWgsToGk);
            meterPoints.add(new TrackPoint(tp.getTime(), c.x, c.y));
        }

        // 按距离断开窗口
        List<List<TrackPoint>> segments = new ArrayList<>();
        List<TrackPoint> currentSegment = new ArrayList<>();
        currentSegment.add(meterPoints.get(0));

        for (int i = 1; i < meterPoints.size(); i++) {
            TrackPoint prev = meterPoints.get(i - 1);
            TrackPoint curr = meterPoints.get(i);

            // 使用米制坐标计算欧几里得距离
            double distance = Math
                    .sqrt(Math.pow(curr.getLon() - prev.getLon(), 2) + Math.pow(curr.getLat() - prev.getLat(), 2));

            if (distance > maxDistanceM) {
                // 距离超过阈值，结束当前段，开始新段
                if (currentSegment.size() >= 2) {
                    segments.add(new ArrayList<>(currentSegment));
                }
                currentSegment.clear();
            }

            currentSegment.add(curr);
        }

        // 添加最后一段
        if (currentSegment.size() >= 2) {
            segments.add(currentSegment);
        }

        if (segments.isEmpty()) {
            return config.geometryFactory.createPolygon();
        }

        // 为每个段构建线缓冲
        List<Geometry> segmentBuffers = new ArrayList<>();
        for (List<TrackPoint> segment : segments) {
            if (segment.size() >= 2) {
                Geometry buffer = buildLineBufferForSegment(segment, widthM);
                if (buffer != null && !buffer.isEmpty()) {
                    segmentBuffers.add(buffer);
                }
            }
        }

        if (segmentBuffers.isEmpty()) {
            return config.geometryFactory.createPolygon();
        }

        // 合并所有缓冲区
        Geometry union = segmentBuffers.get(0);
        for (int i = 1; i < segmentBuffers.size(); i++) {
            union = union.union(segmentBuffers.get(i));
        }

        // 直接返回合并后的缓冲区，不进行凹包平滑
        double originLon = points.get(0).getLon();

        // 转回WGS84
        MathTransform txGkToWgs = getTxGaussToWgsByLon(originLon);
        return JTS.transform(union, txGkToWgs);
    }

    /**
     * 为单个段构建线缓冲
     */
    private Geometry buildLineBufferForSegment(List<TrackPoint> meterSegment, double widthM) {
        if (meterSegment.size() < 2) {
            return null;
        }

        Coordinate[] coords = new Coordinate[meterSegment.size()];
        for (int i = 0; i < meterSegment.size(); i++) {
            TrackPoint tp = meterSegment.get(i);
            coords[i] = new Coordinate(tp.getLon(), tp.getLat());
        }

        LineString line = config.geometryFactory.createLineString(coords);
        Geometry buffer = line.buffer(widthM, config.DEFAULT_BUFFER_QUADRANT);

        return buffer;
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
        double widthM = totalWidthM / 2.0;
        int maxSegmentsToUse = (maxSegments == null || maxSegments <= 0) ? config.DEFAULT_MAX_OUTLINE_SEGMENTS
                : Math.min(maxSegments.intValue(), config.DEFAULT_MAX_OUTLINE_SEGMENTS);
        log.debug("[splitRoad] 开始 参数 原始点数={} 总宽度={}m 返回上限={}", seg.size(), totalWidthM, maxSegmentsToUse);

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

        // 3) 使用距离断开窗口的线缓冲构建轮廓
        log.debug("[splitRoad] 构建线缓冲（距离断开窗口）开始 宽度={}m 输入点数={} 最大段数={}", widthM, workSeg.size(), maxSegmentsToUse);
        long tBuildStart = System.currentTimeMillis();
        Geometry outline = buildOutlineByLineBuffersWithDistanceBreak(workSeg, widthM,
                widthM * config.DISTANCE_BREAK_FACTOR);
        long tBuildEnd = System.currentTimeMillis();
        log.debug("[splitRoad] 构建线缓冲（距离断开窗口）结束 返回 type={} parts={} 耗时={}ms", outline.getGeometryType(),
                (outline instanceof MultiPolygon) ? outline.getNumGeometries() : 1, (tBuildEnd - tBuildStart));

        // 4) 根据maxSegments参数限制返回的多边形数量，并构建parts列表
        log.debug("[splitRoad] 多边形数量限制 开始 最大数量={}", maxSegmentsToUse);
        Geometry finalOutline = outline;
        List<OutlinePart> finalParts = new ArrayList<>();

        if (outline != null && !outline.isEmpty()) {
            List<Geometry> polygons = new ArrayList<>();

            // 收集所有多边形
            if (outline instanceof MultiPolygon) {
                MultiPolygon multiPoly = (MultiPolygon) outline;
                for (int i = 0; i < multiPoly.getNumGeometries(); i++) {
                    polygons.add(multiPoly.getGeometryN(i));
                }
            } else if (outline instanceof Polygon) {
                polygons.add(outline);
            }

            log.debug("[splitRoad] 原始多边形数量={}", polygons.size());

            // 为每个多边形创建OutlinePart（按面积排序）
            List<OutlinePart> allParts = new ArrayList<>();
            for (Geometry polygon : polygons) {
                if (polygon instanceof Polygon) {
                    // 使用正确的calcMu方法计算球面面积
                    double areaMu = calcMu(polygon);
                    String wkt = polygon.toText();

                    // 创建OutlinePart，时间和轨迹点信息暂时留空
                    OutlinePart part = new OutlinePart(
                            polygon,
                            null, // startTime
                            null, // endTime
                            areaMu,
                            wkt);
                    allParts.add(part);
                }
            }

            // 按面积降序排序
            allParts.sort((p1, p2) -> Double.compare(p2.getMu(), p1.getMu()));

            // 选择前maxSegmentsToUse个最大的多边形
            List<OutlinePart> selectedParts = allParts.size() > maxSegmentsToUse
                    ? allParts.subList(0, maxSegmentsToUse)
                    : allParts;

            log.debug("[splitRoad] 保留多边形数量={}", selectedParts.size());

            // 5) 根据最小亩数阈值过滤多边形
            log.debug("[splitRoad] 亩数过滤 开始 最小亩数阈值={}", config.MIN_AREA_MU);
            List<OutlinePart> filteredParts = new ArrayList<>();
            for (OutlinePart part : selectedParts) {
                if (part.getMu() >= config.MIN_AREA_MU) {
                    filteredParts.add(part);
                } else {
                    log.debug("[splitRoad] 过滤掉小面积多边形 面积={}亩 阈值={}亩", part.getMu(), config.MIN_AREA_MU);
                }
            }
            log.debug("[splitRoad] 亩数过滤 结束 过滤前数量={} 过滤后数量={}", selectedParts.size(), filteredParts.size());

            // 重新构建几何结果
            if (filteredParts.size() == 1) {
                finalOutline = filteredParts.get(0).getOutline();
            } else if (filteredParts.size() > 1) {
                GeometryFactory factory = config.geometryFactory;
                Polygon[] selectedPolygons = filteredParts.stream()
                        .map(part -> (Polygon) part.getOutline())
                        .toArray(Polygon[]::new);
                finalOutline = factory.createMultiPolygon(selectedPolygons);
            } else {
                // 如果没有多边形满足亩数要求，返回空多边形
                finalOutline = config.geometryFactory.createPolygon();
            }

            finalParts = filteredParts;

            log.debug("[splitRoad] 多边形数量限制 结束 原始数量={} 保留数量={} 过滤后数量={}",
                    allParts.size(), selectedParts.size(), filteredParts.size());
        }

        // 生成WKT字符串用于输出
        String wkt = finalOutline != null && !finalOutline.isEmpty() ? finalOutline.toText() : null;
        return new SplitRoadResult(finalOutline, finalParts, wkt);
    }

}