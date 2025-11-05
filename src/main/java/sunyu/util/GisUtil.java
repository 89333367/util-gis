package sunyu.util;

import java.time.LocalDateTime;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
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

        // 圆近似细分：让小面积轨迹缓冲区更圆滑（半圆效果）
        private final int DEFAULT_BUFFER_QUADRANT = 2;

        // 作业最大速度阈值（km/h），用于前置速度过滤；可通过 Builder 配置
        private double WORK_MAX_SPEED_KMH = 15.0;
        // 作业最小速度阈值（km/h），用于前置速度过滤；可通过 Builder 配置
        private double MIN_WORK_SPEED_KMH = 1.0;
        // 最小亩数动态阈值（亩），用于 splitRoad 动态过滤小块
        private double MIN_MU_DYNAMIC_THRESHOLD_MU = 0.23;

        // 外缘细长裁剪开关（仅 splitRoad 使用）
        private boolean ENABLE_OUTER_THIN_TRIM = true;
        // 拓扑简化容差系数（外缘细长裁剪用）：控制拓扑保留简化的精度，值越小精度越高
        private double TOPOLOGY_SIMPLIFY_TOLERANCE_FACTOR = 0.2;
        // 拓扑简化容差上限（米）：防止大半径时容差过大
        private double TOPOLOGY_SIMPLIFY_MAX_TOLERANCE_M = 1.0;
        // 外缘细长裁剪半径系数（相对单侧宽度）
        private double THIN_TRIM_RADIUS_FACTOR = 1.5;
        // 宽幅到半径系数的映射表（使用TreeMap支持范围查找）
        private final TreeMap<Double, Double> WIDTH_TO_RADIUS_FACTOR_MAP = new TreeMap<Double, Double>() {
            {
                // 初始化默认映射关系
                put(1.0, 0.45);
                put(1.75, 1.4);
                put(2.5, 1.6);
                put(2.6, 1.8);
                put(2.7, 2.0);
                put(2.8, 2.4);
                put(3.0, 3.0);
            }
        };

        // 线段断裂控制开关（buildSmartLineBuffer 使用）
        private boolean ENABLE_LINE_BREAK = true;
        // 线段断裂距离系数（倍数*单侧宽度），超过则切分会话
        private double LINE_BREAK_FACTOR = 4;

        // 线简化控制开关（buildSmartLineBuffer 使用）
        private boolean ENABLE_LINE_SIMPLIFY = true;
        // 线简化公差系数（倍数*单侧宽度），用于Douglas-Peucker
        private double LINE_SIMPLIFY_TOL_FACTOR = 2;

        // 线缓冲样式（更圆滑边界）：拐角样式、端头样式与mitre限制
        private int BUFFER_JOIN_STYLE = BufferParameters.JOIN_ROUND;
        // 线缓冲的端头样式（圆端），影响线段端点的缓冲形状
        private int BUFFER_END_CAP_STYLE = BufferParameters.CAP_ROUND;
        // 线缓冲的mitre限制值，控制锐角的处理方式，值越小角越圆滑，值越大允许越尖的角
        private double BUFFER_MITRE_LIMIT = 2.0;

        // 兜底方案控制：当splitRoad结果为空或面积过小时启用getOutline重算
        private boolean ENABLE_FALLBACK_TO_OUTLINE = true;
        // 兜底面积阈值（亩）：小于该值时触发兜底方案，使用getOutline重算
        private double FALLBACK_MU_THRESHOLD = 1.5;
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
     * 根据作业宽幅获取对应的外缘细长裁剪半径系数。
     * 
     * <p>
     * 使用TreeMap的floorEntry方法查找小于等于指定宽度的最大键值对，
     * 实现基于宽幅范围的半径系数映射。该系数用于控制外缘细长裁剪的
     * 效果强度，值越大裁剪范围越广。
     * </p>
     * 
     * <p>
     * 映射规则（宽度→半径系数）：
     * <ul>
     * <li>≤0.0米：0.5</li>
     * <li>≤2.5米：1.6</li>
     * <li>≤2.6米：1.8</li>
     * <li>≤2.7米：2.0</li>
     * <li>≤2.8米：2.4</li>
     * <li>≤3.0米：3.0</li>
     * <li>>3.0米：1.5（默认值）</li>
     * </ul>
     * </p>
     * 
     * @param width 作业宽幅（米），必须≥0
     * @return 对应的半径系数，值域[0.5, 3.0]
     * @throws IllegalArgumentException 如果width为负数
     */
    private double getRadiusFactorByWidth(double width) {
        // 参数校验：宽幅不能为负数
        if (width < 0) {
            throw new IllegalArgumentException("作业宽幅不能为负数: " + width);
        }

        // 使用TreeMap.floorEntry查找小于等于width的最大键值对
        // 例如：width=2.55，返回键为2.5的映射（系数1.6）
        Map.Entry<Double, Double> entry = config.WIDTH_TO_RADIUS_FACTOR_MAP.floorEntry(width);
        if (entry != null) {
            return entry.getValue();
        }

        // 当width大于映射表中所有键时，使用默认值1.5
        // 例如：width=5.0，返回THIN_TRIM_RADIUS_FACTOR（1.5）
        return config.THIN_TRIM_RADIUS_FACTOR;
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
     * 构建轨迹线缓冲轮廓 - 支持智能分段的高级版本。
     * 
     * <p>
     * 核心功能：将GPS轨迹点序列转换为指定宽度的面状轮廓，支持基于距离阈值的智能分段。
     * 相比简化版本 {@link #buildSimpleLineBuffer}，本方法提供了：
     * <ul>
     * <li>智能分段：基于距离阈值自动分割轨迹</li>
     * <li>性能优化：支持道格拉斯-普克简化算法</li>
     * <li>详细监控：完整的性能指标记录</li>
     * </ul>
     * </p>
     * 
     * <p>
     * <b>处理流程：</b>
     * </p>
     * <ol>
     * <li>坐标投影：WGS84 → 高斯-克吕格（米制坐标）</li>
     * <li>轨迹分段：基于距离阈值分割轨迹（仅当enableLineBreak=true时）</li>
     * <li>线简化：道格拉斯-普克算法减少点密度（可选）</li>
     * <li>缓冲计算：对每段线段进行指定宽度的缓冲</li>
     * <li>坐标回转：高斯-克吕格 → WGS84</li>
     * </ol>
     * 
     * @param seg             轨迹点列表（WGS84坐标系，建议按时间升序排列）
     * @param widthM          缓冲宽度（米），必须为正数
     * @param enableLineBreak 是否启用基于距离的智能分段（true=启用距离分段，false=整体缓冲）
     * @return 缓冲后的轮廓几何（WGS84坐标系，Polygon或MultiPolygon）
     * @throws IllegalArgumentException 如果参数无效
     * @throws Exception                投影转换或几何运算失败时抛出
     * 
     * @see #buildSimpleLineBuffer(List, double)
     */
    private Geometry buildSmartLineBuffer(List<TrackPoint> seg, double widthM, boolean enableLineBreak)
            throws Exception {
        // 参数校验
        if (seg == null) {
            throw new IllegalArgumentException("轨迹点列表不能为null");
        }
        if (widthM <= 0) {
            throw new IllegalArgumentException("缓冲宽度必须为正数: " + widthM);
        }

        long startTime = System.currentTimeMillis();
        log.trace("[buildSmartLineBuffer] 开始处理 {} 个轨迹点，线缓冲宽度: {} 米", seg.size(), widthM);

        // 快速返回：轨迹点不足无法构成线段
        if (seg.size() < 2) {
            log.debug("[buildSmartLineBuffer] 轨迹点数量不足: {} < 2，返回空几何", seg.size());
            return config.geometryFactory.createGeometryCollection(null);
        }

        // 获取原点经度用于选择合适的高斯-克吕格投影带
        double originLon = seg.get(0).getLon();
        MathTransform txWgsToGk = getTxWgsToGaussByLon(originLon);

        // 计算分段阈值：启用分段时使用配置值，否则禁用分段
        double breakDist = enableLineBreak ? Math.max(widthM, config.LINE_BREAK_FACTOR) : Double.MAX_VALUE;
        log.debug("[buildSmartLineBuffer] 分段控制: enable={}, breakDist={}m (widthM={}m)",
                enableLineBreak, breakDist == Double.MAX_VALUE ? "禁用" : breakDist + "m", widthM);

        // 轨迹分段处理
        List<List<Coordinate>> lines = new ArrayList<>();
        List<Coordinate> currentSegment = new ArrayList<>();
        long convertTime = 0; // 坐标转换耗时统计

        // 遍历轨迹点进行坐标转换和分段
        for (int i = 0; i < seg.size(); i++) {
            TrackPoint p = seg.get(i);

            // WGS84 → 高斯-克吕格坐标转换
            Coordinate src = new Coordinate(p.getLon(), p.getLat());
            Coordinate targetCoord = new Coordinate();
            long t = System.currentTimeMillis();
            JTS.transform(src, targetCoord, txWgsToGk);
            convertTime += System.currentTimeMillis() - t;

            // 距离分段逻辑（仅在启用且当前段非空时执行）
            if (!currentSegment.isEmpty() && enableLineBreak) {
                Coordinate prev = currentSegment.get(currentSegment.size() - 1);
                double dist = Math.hypot(targetCoord.x - prev.x, targetCoord.y - prev.y);

                if (dist > breakDist) {
                    log.trace("[buildSmartLineBuffer] 距离分段: {}m > {}m, 当前段 {} 点完成",
                            dist, breakDist, currentSegment.size());

                    if (currentSegment.size() >= 2) {
                        lines.add(new ArrayList<>(currentSegment)); // 保存当前段
                    }
                    currentSegment.clear(); // 开始新段
                }
            }
            currentSegment.add(targetCoord);
        }

        // 处理最后一段
        if (currentSegment.size() >= 2) {
            lines.add(currentSegment);
            log.debug("[buildSmartLineBuffer] 最后一段: {} 点", currentSegment.size());
        }
        log.debug("[buildSmartLineBuffer] 分段完成: 共 {} 段", lines.size());

        // 几何处理阶段：简化 → 缓冲 → 合并
        long simplifyTime = 0, bufferTime = 0;
        List<LineString> simplifiedLines = new ArrayList<>();

        // 配置缓冲参数
        BufferParameters params = new BufferParameters();
        params.setQuadrantSegments(config.DEFAULT_BUFFER_QUADRANT);
        params.setEndCapStyle(config.BUFFER_END_CAP_STYLE);
        params.setJoinStyle(config.BUFFER_JOIN_STYLE);
        params.setMitreLimit(config.BUFFER_MITRE_LIMIT);

        // 计算简化容差：基于宽度动态调整，确保精度
        double simplifyTolerance = config.ENABLE_LINE_SIMPLIFY
                ? Math.max(0.01, widthM * config.LINE_SIMPLIFY_TOL_FACTOR)
                : 0.0;

        // 处理每个线段：创建折线 → 简化 → 收集
        for (List<Coordinate> coords : lines) {
            if (coords.size() < 2)
                continue; // 跳过无效段

            Coordinate[] coordArray = coords.toArray(new Coordinate[0]);
            LineString line = config.geometryFactory.createLineString(coordArray);

            // 道格拉斯-普克简化（可选）
            long ts = System.currentTimeMillis();
            Geometry simplified = config.ENABLE_LINE_SIMPLIFY
                    ? DouglasPeuckerSimplifier.simplify(line, simplifyTolerance)
                    : line;
            simplifyTime += System.currentTimeMillis() - ts;

            // 处理简化结果（支持MultiLineString）
            if (simplified instanceof LineString) {
                simplifiedLines.add((LineString) simplified);
            } else if (simplified instanceof MultiLineString) {
                MultiLineString mls = (MultiLineString) simplified;
                for (int i = 0; i < mls.getNumGeometries(); i++) {
                    simplifiedLines.add((LineString) mls.getGeometryN(i));
                }
            }
        }

        // 合并所有线段并执行缓冲
        LineString[] lineArray = simplifiedLines.toArray(new LineString[0]);
        MultiLineString multiLine = config.geometryFactory.createMultiLineString(lineArray);

        long tb = System.currentTimeMillis();
        Geometry buffered = BufferOp.bufferOp(multiLine, widthM, params);
        bufferTime += System.currentTimeMillis() - tb;

        // 坐标回转：高斯-克吕格 → WGS84
        long backStart = System.currentTimeMillis();
        MathTransform txGkToWgs = getTxGaussToWgsByLon(originLon);
        Geometry result = JTS.transform(buffered, txGkToWgs);
        long backTime = System.currentTimeMillis() - backStart;

        // 性能统计
        long totalTime = System.currentTimeMillis() - startTime;
        log.debug("[buildSmartLineBuffer] 完成: {}ms (坐标转换:{}ms, 简化:{}ms, 缓冲:{}ms, 回转:{}ms)",
                totalTime, convertTime, simplifyTime, bufferTime, backTime);

        return result;
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
     * 按最小亩数过滤区块（仅对 MultiPolygon 生效）。
     * 规则：
     * - `Polygon`：不进行过滤，直接保留；
     * - `MultiPolygon`：逐面计算亩数，保留面积 ≥ `minMu` 的面；
     * 若全部移除则返回空 `MultiPolygon`；仅剩一个面时返回该 `Polygon`。
     *
     * @param geometry 输入几何（`Polygon`/`MultiPolygon`），可为空
     * @param minMu    最小保留亩数阈值
     * @return 过滤后的几何；类型可能为 `Polygon`、`MultiPolygon` 或原值（非面类型）
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
     * 仅在 `splitRoad` 中启用该裁剪功能。
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

            // 可选：对外环集合做拓扑保留简化以减小缓冲成本（容差可配置，随半径调整）
            double tolThin = Math.min(config.TOPOLOGY_SIMPLIFY_MAX_TOLERANCE_M,
                    radiusM * config.TOPOLOGY_SIMPLIFY_TOLERANCE_FACTOR);
            Geometry shellsSimplified = TopologyPreservingSimplifier.simplify(shellsGeomGk, tolThin);
            if (shellsSimplified == null || shellsSimplified.isEmpty())
                shellsSimplified = shellsGeomGk;

            // 使用配置化缓冲参数，统一执行一次开运算（负缓冲再正缓冲）
            BufferParameters params = new BufferParameters(
                    config.DEFAULT_BUFFER_QUADRANT,
                    BufferParameters.CAP_ROUND,
                    BufferParameters.JOIN_ROUND,
                    2.0);
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
     * 构建区块结果列表
     * 根据几何类型（Polygon或MultiPolygon）生成对应的OutlinePart列表
     * 
     * @param geometry    几何体（Polygon或MultiPolygon）
     * @param trackPoints 轨迹点列表
     * @param totalWidthM 总宽度（米）
     * @return OutlinePart列表
     */
    private List<OutlinePart> buildOutlineParts(Geometry geometry, List<TrackPoint> trackPoints, double totalWidthM) {
        List<OutlinePart> parts = new ArrayList<>();

        if (geometry instanceof Polygon) {
            Polygon poly = (Polygon) geometry;
            OutlinePart part = buildSingleOutlinePart(poly, trackPoints, totalWidthM);
            if (part != null) {
                parts.add(part);
            }
        } else if (geometry instanceof MultiPolygon) {
            MultiPolygon mp = (MultiPolygon) geometry;
            for (int i = 0; i < mp.getNumGeometries(); i++) {
                Polygon poly = (Polygon) mp.getGeometryN(i);
                OutlinePart part = buildSingleOutlinePart(poly, trackPoints, totalWidthM);
                if (part != null) {
                    parts.add(part);
                }
            }
        }

        return parts;
    }

    /**
     * 构建单个区块结果
     * 
     * @param polygon     多边形
     * @param trackPoints 轨迹点列表
     * @param totalWidthM 总宽度（米）
     * @return OutlinePart对象，如果构建失败返回null
     */
    private OutlinePart buildSingleOutlinePart(Polygon polygon, List<TrackPoint> trackPoints, double totalWidthM) {
        try {
            // 计算面积
            double mu = calcMu(polygon);

            // 生成WKT
            String wkt = toWkt(polygon);

            // 统计轮廓内的轨迹点
            List<TrackPoint> pointsInPolygon = filterPointsInPolygon(trackPoints, polygon);

            // 提取时间范围
            LocalDateTime startTime = null;
            LocalDateTime endTime = null;
            if (!pointsInPolygon.isEmpty()) {
                startTime = pointsInPolygon.get(0).getTime();
                endTime = pointsInPolygon.get(pointsInPolygon.size() - 1).getTime();
            }

            return new OutlinePart(polygon, startTime, endTime, mu, wkt, pointsInPolygon)
                    .setTotalWidthM(totalWidthM);
        } catch (Exception e) {
            log.warn("构建单个区块结果失败", e);
            return null;
        }
    }

    /**
     * 过滤多边形内的轨迹点
     * 
     * @param trackPoints 轨迹点列表
     * @param polygon     多边形
     * @return 多边形内的轨迹点列表
     */
    private List<TrackPoint> filterPointsInPolygon(List<TrackPoint> trackPoints, Polygon polygon) {
        if (trackPoints == null || trackPoints.isEmpty()) {
            return new ArrayList<>();
        }

        PreparedGeometry preparedPoly = PreparedGeometryFactory.prepare(polygon);
        Envelope env = polygon.getEnvelopeInternal();

        return trackPoints.stream()
                .filter(p -> {
                    Coordinate c = new Coordinate(p.getLon(), p.getLat());
                    // 先进行边界框快速过滤
                    if (!env.contains(c)) {
                        return false;
                    }
                    // 再进行精确包含判断
                    Geometry point = config.geometryFactory.createPoint(c);
                    return preparedPoly.contains(point);
                })
                .collect(Collectors.toList());
    }

    /**
     * 判断是否需要触发兜底方案
     * 当结果为空或总亩数小于阈值时触发兜底
     * 
     * @param result 拆分结果
     * @return 是否需要兜底
     */
    private boolean shouldUseFallback(SplitRoadResult result) {
        // 如果配置关闭了兜底方案，直接返回false
        if (!config.ENABLE_FALLBACK_TO_OUTLINE) {
            return false;
        }

        if (result == null) {
            return true;
        }

        // 结果为空
        if (result.getParts() == null || result.getParts().isEmpty()) {
            return true;
        }

        // 总亩数过小（小于配置的兜底面积阈值）
        double totalMu = result.getMu();
        if (totalMu < config.FALLBACK_MU_THRESHOLD) {
            return true;
        }

        return false;
    }

    /**
     * 应用兜底方案，使用getOutline方法重新计算
     * 
     * @param trackPoints 轨迹点列表
     * @param totalWidthM 总宽度（米）
     * @return 兜底方案结果
     * @throws Exception 计算异常
     */
    private SplitRoadResult applyFallback(List<TrackPoint> trackPoints, double totalWidthM) throws Exception {
        try {
            // 使用getOutline方法重新计算（返回OutlinePart）
            OutlinePart outlinePart = getOutline(trackPoints, totalWidthM);

            // 构建单个区块的列表
            List<OutlinePart> parts = new ArrayList<>();
            if (outlinePart != null) {
                parts.add(outlinePart);
            }

            // 获取几何体和WKT
            Geometry outline = outlinePart.getOutline();
            String outlineWkt = outlinePart.getWkt();

            // 构建结果对象
            SplitRoadResult result = new SplitRoadResult(outline, parts, outlineWkt)
                    .setTotalWidthM(totalWidthM);
            result.calculateTimeRange();

            log.info("[splitRoad] 兜底方案执行成功: 区块数={}, 总亩数={}",
                    result.getParts().size(), result.getMu());

            return result;

        } catch (Exception e) {
            log.trace("[splitRoad] 兜底方案执行失败", e);

            // 兜底方案失败时，返回包含空几何的SplitRoadResult，避免调用者报错
            Geometry emptyGeometry = config.geometryFactory.createGeometryCollection(null);
            String emptyWkt = "GEOMETRYCOLLECTION EMPTY";
            List<OutlinePart> emptyParts = new ArrayList<>();

            SplitRoadResult emptyResult = new SplitRoadResult(emptyGeometry, emptyParts, emptyWkt)
                    .setTotalWidthM(totalWidthM);
            emptyResult.calculateTimeRange();

            log.warn("[splitRoad] 兜底方案失败，返回空结果: 错误={}", e.getMessage());
            return emptyResult;
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
        return new OutlinePart(outline, startTime, endTime, mu, wkt, points)
                .setTotalWidthM(totalWidthM);
    }

    /**
     * 基于轨迹点与总宽度生成轮廓，并按面积保留 Top-N 最大区块。
     * 
     * 处理流程：
     * 1. 异常点过滤与排序
     * 2. 速度范围过滤
     * 3. 线缓冲构建轮廓
     * 4. 外缘细长条裁剪（可选）
     * 5. 保留最大区块并过滤小面积
     * 6. 生成区块结果
     * 
     * @param seg         轨迹点列表（WGS84），至少3个有效点
     * @param totalWidthM 作业总宽度（米），左右合计，必须非负
     * @return 拆分结果，包含几何体和区块信息
     * @throws Exception 坐标转换或几何运算异常
     */
    public SplitRoadResult splitRoad(List<TrackPoint> seg, double totalWidthM) throws Exception {
        return splitRoad(seg, totalWidthM, config.DEFAULT_MAX_OUTLINE_SEGMENTS);
    }

    /**
     * 基于轨迹点与总宽度生成轮廓，并按面积保留不超过 maxSegments 的最大区块。
     * 
     * 处理流程：
     * 1. 异常点过滤与排序
     * 2. 速度范围过滤
     * 3. 线缓冲构建轮廓
     * 4. 外缘细长条裁剪（可选）
     * 5. 保留最大区块并过滤小面积
     * 6. 生成区块结果
     * 
     * @param seg         轨迹点列表（WGS84），至少3个有效点
     * @param totalWidthM 作业总宽度（米），左右合计，必须非负
     * @param maxSegments 返回的区块数量上限；为null或<=0时取默认值
     * @return 拆分结果，包含几何体和区块信息
     * @throws Exception 坐标转换或几何运算异常
     */
    public SplitRoadResult splitRoad(List<TrackPoint> seg, double totalWidthM, Integer maxSegments) throws Exception {
        long startTime = System.currentTimeMillis();
        double halfWidth = totalWidthM / 2.0;
        int segmentLimit = (maxSegments == null || maxSegments <= 0) ? config.DEFAULT_MAX_OUTLINE_SEGMENTS
                : maxSegments;

        log.debug("[splitRoad] 开始处理: 点数={}, 宽度={}m, 上限={}", seg.size(), totalWidthM, segmentLimit);

        // 1) 过滤异常点（越界与坐标为0的点），按时间升序
        long filterStart = System.currentTimeMillis();
        List<TrackPoint> sortedPoints = filterAndSortTrackPoints(seg);
        long filterEnd = System.currentTimeMillis();
        log.debug("[splitRoad] 异常点过滤完成: 有效点数={}, 耗时={}ms", sortedPoints.size(), (filterEnd - filterStart));

        // 2) 速度过滤：删除小于等于1km/h或大于等于15km/h的点
        long speedStart = System.currentTimeMillis();
        List<TrackPoint> validPoints = filterBySpeedRange(sortedPoints, config.MIN_WORK_SPEED_KMH,
                config.WORK_MAX_SPEED_KMH);
        long speedEnd = System.currentTimeMillis();
        int removedCount = sortedPoints.size() - validPoints.size();
        log.debug("[splitRoad] 速度过滤完成: 移除点数={}, 剩余点数={}, 耗时={}ms",
                removedCount, validPoints.size(), (speedEnd - speedStart));

        // 3) 构建线缓冲轮廓
        long buildStart = System.currentTimeMillis();
        Geometry outline = buildSmartLineBuffer(validPoints, halfWidth, config.ENABLE_LINE_BREAK);
        long buildEnd = System.currentTimeMillis();
        int outlineParts = (outline instanceof MultiPolygon) ? outline.getNumGeometries() : 1;
        log.debug("[splitRoad] 线缓冲构建完成: 类型={}, 区块数={}, 耗时={}ms",
                outline.getGeometryType(), outlineParts, (buildEnd - buildStart));

        // 4) 外缘细长条裁剪（可选）
        if (config.ENABLE_OUTER_THIN_TRIM) {
            long thinStart = System.currentTimeMillis();
            double radiusFactor = getRadiusFactorByWidth(totalWidthM);
            double radiusM = Math.max(0.1, radiusFactor);
            int partsBeforeThin = (outline instanceof MultiPolygon) ? outline.getNumGeometries() : 1;

            log.debug("[splitRoad] 外缘细长裁剪: 半径={}m (系数={}), 输入区块数={}", radiusM, radiusFactor, partsBeforeThin);
            Geometry outlineAfterThin = trimOuterThinStrips(outline, radiusM);

            int partsAfterThin = (outlineAfterThin instanceof MultiPolygon) ? outlineAfterThin.getNumGeometries() : 1;
            int removedThin = Math.max(0, partsBeforeThin - partsAfterThin);
            outline = outlineAfterThin;

            long thinEnd = System.currentTimeMillis();
            log.debug("[splitRoad] 外缘细长裁剪完成: 移除区块数={}, 剩余区块数={}, 耗时={}ms",
                    removedThin, partsAfterThin, (thinEnd - thinStart));
        }

        // 5) 保留最大区块并过滤小面积
        log.debug("[splitRoad] 区块筛选: 上限={}, 输入类型={}, 区块数={}", segmentLimit, outline.getGeometryType(), outlineParts);

        long keepStart = System.currentTimeMillis();
        Geometry kept = keepLargestPolygons(outline, segmentLimit);
        long keepEnd = System.currentTimeMillis();
        int partsAfterKeep = (kept instanceof MultiPolygon) ? kept.getNumGeometries() : 1;
        log.debug("[splitRoad] 保留大区块完成: 输出区块数={}, 耗时={}ms", partsAfterKeep, (keepEnd - keepStart));

        long areaStart = System.currentTimeMillis();
        Geometry trimmed = removeSmallMuPolygons(kept, config.MIN_MU_DYNAMIC_THRESHOLD_MU);
        long areaEnd = System.currentTimeMillis();
        int partsAfterArea = (trimmed instanceof MultiPolygon) ? trimmed.getNumGeometries() : 1;
        int removedByArea = partsAfterKeep - partsAfterArea;
        log.debug("[splitRoad] 小面积过滤完成: 移除区块数={}, 剩余区块数={}, 耗时={}ms",
                removedByArea, partsAfterArea, (areaEnd - areaStart));

        // 6) 生成区块结果
        List<OutlinePart> parts = buildOutlineParts(trimmed, validPoints, totalWidthM);

        long wktStart = System.currentTimeMillis();
        String outlineWkt = toWkt(trimmed);
        long wktEnd = System.currentTimeMillis();
        log.trace("[splitRoad] WKT生成完成: 长度={}, 耗时={}ms", outlineWkt.length(), (wktEnd - wktStart));

        long totalEnd = System.currentTimeMillis();
        log.debug("[splitRoad] 处理完成: 区块数={}, 总耗时={}ms", parts.size(), (totalEnd - startTime));

        SplitRoadResult result = new SplitRoadResult(trimmed, parts, outlineWkt).setTotalWidthM(totalWidthM);
        result.calculateTimeRange();

        // 兜底方案：结果为空或亩数过小时使用getOutline重算
        if (shouldUseFallback(result)) {
            log.warn("[splitRoad] 触发兜底方案: 区块数={}, 总亩数={}",
                    result.getParts().size(), result.getMu());
            return applyFallback(validPoints, totalWidthM);
        }

        return result;
    }

}