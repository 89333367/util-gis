package sunyu.util;

import cn.hutool.core.bean.BeanUtil;
import cn.hutool.log.Log;
import cn.hutool.log.LogFactory;
import org.geotools.api.referencing.crs.CoordinateReferenceSystem;
import org.geotools.api.referencing.operation.MathTransform;
import org.geotools.geojson.geom.GeometryJSON;
import org.geotools.geometry.jts.JTS;
import org.geotools.referencing.CRS;
import org.locationtech.jts.algorithm.hull.ConcaveHull;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.LineString;
import org.locationtech.jts.simplify.DouglasPeuckerSimplifier;
import sunyu.util.pojo.CoordinatePoint;
import sunyu.util.pojo.TrackPoint;

import java.util.Comparator;
import java.util.List;
import java.util.concurrent.ConcurrentHashMap;
import java.util.stream.Collectors;

/**
 * GIS工具类，用于轨迹处理、空间计算等
 *
 * @author SunYu
 */
public class GisUtil implements AutoCloseable {
    // 日志记录器，用于记录工具类的运行状态和调试信息
    private final Log log = LogFactory.get();
    // 配置参数，包含各种常量和默认值
    private final Config config;

    /**
     * 获取构建器实例，使用构建器模式创建GisUtil对象
     * 构建器模式允许逐步构建复杂对象，提高代码可读性和灵活性
     *
     * @return Builder对象，用于构建GisUtil实例
     */
    public static Builder builder() {
        return new Builder();
    }

    /**
     * 私有构造函数，通过Builder创建GisUtil实例
     * 采用私有构造函数防止直接实例化，强制使用Builder模式
     *
     * @param config 配置参数，包含各种常量和默认值
     */
    private GisUtil(Config config) {
        // 记录工具类构建开始日志
        log.info("[构建工具类] 开始");
        // 其他初始化语句（预留扩展点）
        // 记录工具类构建结束日志
        log.info("[构建工具类] 结束");
        // 保存配置参数引用
        this.config = config;
    }

    /**
     * 内部配置类，定义一些常量
     * 使用内部类封装配置参数，提高代码组织性和封装性
     */
    private static class Config {
        // 每平方米对应的mu单位（面积单位），用于将平方米转换为亩
        // 1亩 = 666.666...平方米，所以 1平方米 = 1/666.666... ≈ 0.0015000015
        private final double MU_PER_SQ_METER = 0.0015000015;

        // WGS84坐标系的EPSG代码，用于定义地理坐标系统
        private final String WGS84 = "EPSG:4326";

        // 地球半径（米），用于Haversine公式计算两点间距离
        private final double R = 6371000;

        // 长半轴（椭球体的赤道半径）
        private final double semiMajorAxis = 6378245.0;

        // 椭球体偏心率的平方
        private final double eccentricitySquared = 0.00669342162296594323;

        // CRS缓存，避免重复解析WKT
        private final ConcurrentHashMap<String, CoordinateReferenceSystem> crsCache = new ConcurrentHashMap<>();

        // GeometryFactory缓存，避免重复创建
        private final GeometryFactory geometryFactory = new GeometryFactory();
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
        // 清理缓存
        config.crsCache.clear();
        // 回收各种资源（预留扩展点）
    }

    /**
     * 从缓存中获取或创建CoordinateReferenceSystem
     *
     * @param wkt WKT字符串
     *
     * @return CoordinateReferenceSystem对象
     *
     */
    private CoordinateReferenceSystem getCachedCRS(String wkt) {
        return config.crsCache.computeIfAbsent(wkt, key -> {
            try {
                return CRS.parseWKT(key);
            } catch (Exception e) {
                throw new RuntimeException("Failed to parse CRS WKT", e);
            }
        });
    }


    /**
     * 根据经度选择合适的坐标参考系统（CRS）
     * 不同区域使用不同的投影系统以提高计算精度
     *
     * @param lon 经度值，用于判断所在区域
     *
     * @return CoordinateReferenceSystem对象，对应区域的最佳投影系统
     *
     */
    private CoordinateReferenceSystem pickCrs(double lon) {
        // 中国区域使用CGCS2000 3度分带投影（经度72°至138°）
        // CGCS2000是中国2000国家大地坐标系，3度分带适用于中国大部分地区
        if (lon >= 72 && lon <= 138) {
            return getCachedCRS(
                    "PROJCS[\"CGCS2000 / 3-degree Gauss-Kruger zone 39N\","
                            + "GEOGCS[\"CGCS2000\","
                            + "DATUM[\"CGCS2000\",SPHEROID[\"CGCS2000\",6378137,298.257222101]],"
                            + "PRIMEM[\"Greenwich\",0],UNIT[\"degree\",0.0174532925199433]],"
                            + "PROJECTION[\"Transverse_Mercator\"],"
                            + "PARAMETER[\"central_meridian\",117],"
                            + "PARAMETER[\"scale_factor\",1],"
                            + "PARAMETER[\"false_easting\",39500000],"
                            + "PARAMETER[\"false_northing\",0],"
                            + "UNIT[\"metre\",1]]");
        }
        // 北美区域使用NAD83 Albers投影（经度-140°至-50°）
        // Albers投影是一种等积圆锥投影，适用于北美大陆
        if (lon >= -140 && lon <= -50) {
            return getCachedCRS(
                    "PROJCS[\"NAD83 / Conus Albers\","
                            + "GEOGCS[\"NAD83\","
                            + "DATUM[\"NAD83\",SPHEROID[\"GRS 1980\",6378137,298.257222101]],"
                            + "PRIMEM[\"Greenwich\",0],UNIT[\"degree\",0.0174532925199433]],"
                            + "PROJECTION[\"Albers_Conic_Equal_Area\"],"
                            + "PARAMETER[\"standard_parallel_1\",29.5],"
                            + "PARAMETER[\"standard_parallel_2\",45.5],"
                            + "PARAMETER[\"central_meridian\",-96],"
                            + "PARAMETER[\"latitude_of_origin\",37.5],"
                            + "UNIT[\"metre\",1]]");
        }
        // 欧洲区域使用ETRS89 LAEA投影（经度-30°至60°）
        // LAEA（Lambert Azimuthal Equal-Area）投影是一种等积方位投影
        if (lon >= -30 && lon <= 60) {
            return getCachedCRS(
                    "PROJCS[\"ETRS89-extended / LAEA Europe\","
                            + "GEOGCS[\"ETRS89\","
                            + "DATUM[\"ETRS89\",SPHEROID[\"GRS 1980\",6378137,298.257222101]],"
                            + "PRIMEM[\"Greenwich\",0],UNIT[\"degree\",0.0174532925199433]],"
                            + "PROJECTION[\"Lambert_Azimuthal_Equal_Area\"],"
                            + "PARAMETER[\"latitude_of_origin\",52],"
                            + "PARAMETER[\"central_meridian\",10],"
                            + "UNIT[\"metre\",1]]");
        }
        // 其他区域使用世界Mollweide投影（全球等积伪圆柱投影）
        // Mollweide投影是一种等积投影，适用于全球范围的地图显示
        return getCachedCRS(
                "PROJCS[\"World_Mollweide\","
                        + "GEOGCS[\"WGS 84\","
                        + "DATUM[\"WGS_1984\",SPHEROID[\"WGS 84\",6378137,298.257223563]],"
                        + "PRIMEM[\"Greenwich\",0],UNIT[\"degree\",0.0174532925199433]],"
                        + "PROJECTION[\"Mollweide\"],PARAMETER[\"central_meridian\",0],"
                        + "UNIT[\"metre\",1]]");
    }

    /**
     * 构建轨迹轮廓
     * 通过轨迹点生成带宽度的轮廓多边形，用于面积计算等操作
     *
     * @param seg            轨迹点列表，至少需要3个点
     * @param leftWidthM     左侧宽度（米），轨迹线左侧的扩展距离
     * @param rightWidthM    右侧宽度（米），轨迹线右侧的扩展距离
     * @param simplifyM      简化阈值（米），用于轨迹点简化以提高性能
     *                       推荐值范围：
     *                       低精度场景：10-50米（如农田作业粗略估算）
     *                       中等精度场景：5-10米（一般农业机械作业）
     *                       高精度场景：1-5米（精细化农业作业）
     * @param maxEdgeLengthM 最大边缘长度（米），控制轮廓的精细程度
     *                       推荐值范围：
     *                       粗略轮廓：50-100米
     *                       标准轮廓：20-50米
     *                       精细轮廓：5-20米
     *
     * @return Geometry对象，表示生成的轮廓多边形
     *
     * @throws Exception 可能抛出异常，如坐标转换错误或几何操作失败
     */
    private Geometry buildOutline(List<TrackPoint> seg,
                                  double leftWidthM,
                                  double rightWidthM,
                                  double simplifyM,
                                  double maxEdgeLengthM) throws Exception {
        // 如果轨迹点少于3个，无法构成有效几何形状，抛出异常
        // 添加输入验证
        if (seg == null || seg.size() < 3) {
            throw new IllegalArgumentException("轨迹段至少需要3个点");
        }
        if (leftWidthM < 0 || rightWidthM < 0 || simplifyM < 0 || maxEdgeLengthM < 0) {
            throw new IllegalArgumentException("所有参数必须为非负数");
        }

        // 先按时间戳排序，确保轨迹点的时序正确
        List<TrackPoint> sortedSeg = seg.stream()
                .sorted(Comparator.comparing(TrackPoint::getTime))
                .collect(Collectors.toList());
        // 将排序后的轨迹点转换为坐标数组
        Coordinate[] coords = sortedSeg.stream()
                .map(p -> new Coordinate(p.getLon(), p.getLat()))
                .toArray(Coordinate[]::new);
        // 创建原始线段，由所有轨迹点连接而成
        LineString rawLine = config.geometryFactory.createLineString(coords);
        // 使用Douglas-Peucker算法简化线段，减少点数提高性能
        // 使用起点纬度计算更精确的简化阈值
        double degPerMeter = degreesPerMeterAtLat(seg.get(0).getLat());
        LineString simple = (LineString) DouglasPeuckerSimplifier.simplify(rawLine, simplifyM * degPerMeter);

        // 定义源坐标系统（WGS84地理坐标系）
        CoordinateReferenceSystem src = CRS.decode(config.WGS84, true);
        // 根据第一个点的经度选择最佳目标坐标系统（投影坐标系）
        CoordinateReferenceSystem tgt = pickCrs(seg.get(0).getLon()); // 使用原始经度值而不是投影后的x值
        // 获取从源坐标系到目标坐标系的数学变换对象
        MathTransform tx = CRS.findMathTransform(src, tgt, true);
        // 对简化后的线段进行坐标转换，从地理坐标转为投影坐标
        LineString proj = (LineString) JTS.transform(simple, tx);

        // 创建单侧缓冲区然后合并，而不是创建两个独立的缓冲区
        // 这样可以提高性能并减少几何错误
        Geometry strip;
        // 直接使用 equals 判定左右宽度相等（如果业务允许精确比较）
        if (leftWidthM == rightWidthM) {
            strip = proj.buffer(leftWidthM);
        } else {
            // 优化：先创建联合缓冲区，减少几何操作次数
            Geometry leftBuffer = proj.buffer(leftWidthM);
            Geometry rightBuffer = proj.buffer(rightWidthM);
            strip = leftBuffer.union(rightBuffer);
        }

        // 使用凹包算法生成轮廓，maxEdgeLenM控制轮廓边缘的最大长度
        return ConcaveHull.concaveHullByLength(strip, maxEdgeLengthM);
    }


    /**
     * 计算轨迹段的面积（mu单位）
     * 通过构建轨迹轮廓并计算面积，结果以亩为单位
     *
     * @param seg            轨迹段，至少需要3个点
     * @param leftWidthM     左侧宽度（米），轨迹线左侧的扩展距离
     * @param rightWidthM    右侧宽度（米），轨迹线右侧的扩展距离
     * @param simplifyM      简化阈值（米），用于轨迹点简化以提高性能
     *                       推荐值范围：
     *                       低精度场景：10-50米（如农田作业粗略估算）
     *                       中等精度场景：5-10米（一般农业机械作业）
     *                       高精度场景：1-5米（精细化农业作业）
     * @param maxEdgeLengthM 最大边缘长度（米），控制轮廓的精细程度
     *                       推荐值范围：
     *                       粗略轮廓：50-100米
     *                       标准轮廓：20-50米
     *                       精细轮廓：5-20米
     *
     * @return 面积（mu），以亩为单位，保留3位小数
     *
     * @throws Exception 可能抛出异常，如坐标转换错误或几何操作失败
     */
    public double calcMu(List<TrackPoint> seg,
                         double leftWidthM,
                         double rightWidthM,
                         double simplifyM,
                         double maxEdgeLengthM) throws Exception {
        // 构建轨迹轮廓几何对象
        Geometry outline = buildOutline(seg, leftWidthM, rightWidthM, simplifyM, maxEdgeLengthM);
        return calcMu(outline);
    }


    /**
     * 计算轮廓面积(mu单位)
     *
     * @param outline 轮廓
     *
     * @return 面积（mu），以亩为单位，保留3位小数
     */
    public double calcMu(Geometry outline) {
        // 计算面积并转换为mu单位（亩）
        // g.getArea() 获取几何对象面积（平方米）
        // * config.MU_PER_SQ_METER 转换为亩
        // * 1000.0 和 / 1000.0 实现保留3位小数的四舍五入
        return Math.round(outline.getArea() * config.MU_PER_SQ_METER * 1000.0) / 1000.0;
    }


    /**
     * 将轨迹段转换为WKT格式
     * 生成轨迹轮廓的WKT（Well-Known Text）表示，便于GIS软件处理
     *
     * @param seg            轨迹段，至少需要3个点
     * @param leftWidthM     左侧宽度（米），轨迹线左侧的扩展距离
     * @param rightWidthM    右侧宽度（米），轨迹线右侧的扩展距离
     * @param simplifyM      简化阈值（米），用于轨迹点简化以提高性能
     *                       推荐值范围：
     *                       低精度场景：10-50米（如农田作业粗略估算）
     *                       中等精度场景：5-10米（一般农业机械作业）
     *                       高精度场景：1-5米（精细化农业作业）
     * @param maxEdgeLengthM 最大边缘长度（米），控制轮廓的精细程度
     *                       推荐值范围：
     *                       粗略轮廓：50-100米
     *                       标准轮廓：20-50米
     *                       精细轮廓：5-20米
     *
     * @return WKT字符串，表示轨迹轮廓的几何形状
     *
     * @throws Exception 可能抛出异常，如坐标转换错误或几何操作失败
     */
    public String toWkt(List<TrackPoint> seg,
                        double leftWidthM,
                        double rightWidthM,
                        double simplifyM,
                        double maxEdgeLengthM) throws Exception {
        // 构建轨迹轮廓几何对象
        Geometry g = buildOutline(seg, leftWidthM, rightWidthM, simplifyM, maxEdgeLengthM);
        // 转换回WGS84坐标系统，确保输出为标准地理坐标
        // 使用第一个轨迹点的原始经度来选择正确的投影坐标系
        CoordinateReferenceSystem src = pickCrs(seg.get(0).getLon());  // 使用原始经度值
        CoordinateReferenceSystem tgt = CRS.decode(config.WGS84, true);
        MathTransform tx = CRS.findMathTransform(src, tgt, true);
        Geometry wgs = JTS.transform(g, tx);
        // 返回WKT字符串表示
        return wgs.toText();
    }

    /**
     * 使用haversine公式计算两点间距离
     * haversine公式用于计算球面上两点间的最短距离（大圆距离）
     *
     * @param p1 第一个点，包含经纬度信息
     * @param p2 第二个点，包含经纬度信息
     *
     * @return 距离（米），两点间的地理距离
     */
    private double haversine(CoordinatePoint p1, CoordinatePoint p2) {
        // 计算纬度差并转换为弧度
        double dLat = Math.toRadians(p2.getLat() - p1.getLat());
        // 计算经度差并转换为弧度
        double dLon = Math.toRadians(p2.getLon() - p1.getLon());
        // haversine公式的核心计算部分
        // a = sin²(Δlat/2) + cos(lat1) * cos(lat2) * sin²(Δlon/2)
        double a = Math.sin(dLat / 2) * Math.sin(dLat / 2) +
                Math.cos(Math.toRadians(p1.getLat())) *           // 第一点纬度的余弦值
                        Math.cos(Math.toRadians(p2.getLat())) *   // 第二点纬度的余弦值
                        Math.sin(dLon / 2) * Math.sin(dLon / 2);  // 经度差的正弦平方
        // 计算大圆距离，公式为：距离 = 2 * R * arcsin(√a)
        // 其中R为地球半径，a为上面计算的值
        return config.R * 2 * Math.atan2(Math.sqrt(a), Math.sqrt(1 - a));
    }

    /**
     * 计算加密参数
     *
     * @param lon 经度
     * @param lat 纬度
     *
     * @return 加密偏移量数组，[0]为经度偏移，[1]为纬度偏移
     */
    private double[] delta(double lon, double lat) {
        double dLat = transformLat(lon - 105.0, lat - 35.0);
        double dLon = transformLon(lon - 105.0, lat - 35.0);
        double radLat = lat / 180.0 * Math.PI;
        double magic = Math.sin(radLat);
        magic = 1 - config.eccentricitySquared * magic * magic;
        double sqrtMagic = Math.sqrt(magic);
        // 添加保护条件
        if (magic * sqrtMagic < 1e-10) {
            throw new IllegalArgumentException("纬度值导致除零错误");
        }
        dLat = (dLat * 180.0) / ((config.semiMajorAxis * (1 - config.eccentricitySquared)) / (magic * sqrtMagic) * Math.PI);
        dLon = (dLon * 180.0) / (config.semiMajorAxis / sqrtMagic * Math.cos(radLat) * Math.PI);
        return new double[]{dLon, dLat};
    }


    /**
     * 转换纬度
     *
     * @param x 经度差值
     * @param y 纬度差值
     *
     * @return 转换后的纬度值
     */
    private double transformLat(double x, double y) {
        double ret = -100.0 + 2.0 * x + 3.0 * y + 0.2 * y * y + 0.1 * x * y + 0.2 * Math.sqrt(Math.abs(x));
        ret += (20.0 * Math.sin(6.0 * x * Math.PI) + 20.0 * Math.sin(2.0 * x * Math.PI)) * 2.0 / 3.0;
        ret += (20.0 * Math.sin(y * Math.PI) + 40.0 * Math.sin(y / 3.0 * Math.PI)) * 2.0 / 3.0;
        ret += (160.0 * Math.sin(y / 12.0 * Math.PI) + 320 * Math.sin(y * Math.PI / 30.0)) * 2.0 / 3.0;
        return ret;
    }

    /**
     * 转换经度
     *
     * @param x 经度差值
     * @param y 纬度差值
     *
     * @return 转换后的经度值
     */
    private double transformLon(double x, double y) {
        double ret = 300.0 + x + 2.0 * y + 0.1 * x * x + 0.1 * x * y + 0.1 * Math.sqrt(Math.abs(x));
        ret += (20.0 * Math.sin(6.0 * x * Math.PI) + 20.0 * Math.sin(2.0 * x * Math.PI)) * 2.0 / 3.0;
        ret += (20.0 * Math.sin(x * Math.PI) + 40.0 * Math.sin(x / 3.0 * Math.PI)) * 2.0 / 3.0;
        ret += (150.0 * Math.sin(x / 12.0 * Math.PI) + 300.0 * Math.sin(x / 30.0 * Math.PI)) * 2.0 / 3.0;
        return ret;
    }

    /**
     * 根据纬度计算每米对应的经度变化（度）
     * 地球是椭球体，不同纬度上每米对应的经度变化不同
     *
     * @param lat 纬度
     *
     * @return 每米对应的经度变化（度）
     */
    private double degreesPerMeterAtLat(double lat) {
        // 在给定纬度上，地球的半径会变小
        double latRad = Math.toRadians(lat);
        // 使用简单的球面近似计算
        double metersPerDegree = 2 * Math.PI * config.R * Math.cos(latRad) / 360;

        // 防止除零错误，当cos(latRad)接近0时（接近两极）
        if (Math.abs(metersPerDegree) < 1e-10) {
            return 1.0 / (2 * Math.PI * config.R / 360);
        }

        return 1.0 / metersPerDegree;
    }


    /**
     * 判断两个几何形状是否拓扑相等
     *
     * @param g1 第一个几何形状
     * @param g2 第二个几何形状
     *
     * @return 如果两个几何形状拓扑相等则返回true，否则返回false
     */
    public boolean equals(Geometry g1, Geometry g2) {
        return g1.equals(g2);
    }

    /**
     * 判断两个几何形状是否没有共有点（脱节）
     *
     * @param g1 第一个几何形状
     * @param g2 第二个几何形状
     *
     * @return 如果两个几何形状没有共有点则返回true，否则返回false
     */
    public boolean disjoint(Geometry g1, Geometry g2) {
        return g1.disjoint(g2);
    }

    /**
     * 判断两个几何形状是否至少有一个共有点（相交）
     *
     * @param g1 第一个几何形状
     * @param g2 第二个几何形状
     *
     * @return 如果两个几何形状至少有一个共有点则返回true，否则返回false
     */
    public boolean intersects(Geometry g1, Geometry g2) {
        return g1.intersects(g2);
    }

    /**
     * 判断两个几何形状是否有至少一个公共的边界点，但是没有内部点（接触）
     *
     * @param g1 第一个几何形状
     * @param g2 第二个几何形状
     *
     * @return 如果两个几何形状接触则返回true，否则返回false
     */
    public boolean touches(Geometry g1, Geometry g2) {
        return g1.touches(g2);
    }

    /**
     * 判断两个几何形状是否共享一些但不是所有的内部点（交叉）
     *
     * @param g1 第一个几何形状
     * @param g2 第二个几何形状
     *
     * @return 如果两个几何形状交叉则返回true，否则返回false
     */
    public boolean crosses(Geometry g1, Geometry g2) {
        return g1.crosses(g2);
    }

    /**
     * 判断几何形状A是否完全在几何形状B内部（内含）
     *
     * @param g1 几何形状A
     * @param g2 几何形状B
     *
     * @return 如果几何形状A在几何形状B内部则返回true，否则返回false
     */
    public boolean within(Geometry g1, Geometry g2) {
        return g1.within(g2);
    }

    /**
     * 判断几何形状A是否包含几何形状B（包含）
     *
     * @param g1 几何形状A
     * @param g2 几何形状B
     *
     * @return 如果几何形状A包含几何形状B则返回true，否则返回false
     */
    public boolean contains(Geometry g1, Geometry g2) {
        return g1.contains(g2);
    }

    /**
     * 判断两个几何形状是否共享一部分但不是所有的公共点，而且相交处有他们自己相同的区域（重叠）
     *
     * @param g1 第一个几何形状
     * @param g2 第二个几何形状
     *
     * @return 如果两个几何形状重叠则返回true，否则返回false
     */
    public boolean overlaps(Geometry g1, Geometry g2) {
        return g1.overlaps(g2);
    }

    /**
     * 判断点是否在矩形内
     * <p>
     * 注意：在地理坐标系中，经度表示东西方向（西经为负，东经为正），纬度表示南北方向（南纬为负，北纬为正）
     * 所以矩形的"左下角"是(minLon, minLat)，"右上角"是(maxLon, maxLat)
     *
     * @param pointLon   点的经度（X坐标）
     * @param pointLat   点的纬度（Y坐标）
     * @param rectMinLon 矩形最小经度（西边界），对应矩形左边界
     * @param rectMinLat 矩形最小纬度（南边界），对应矩形下边界
     * @param rectMaxLon 矩形最大经度（东边界），对应矩形右边界
     * @param rectMaxLat 矩形最大纬度（北边界），对应矩形上边界
     *
     * @return 如果点在矩形内则返回true，否则返回false
     */
    public boolean inRectangle(double pointLon, double pointLat,
                               double rectMinLon, double rectMinLat,
                               double rectMaxLon, double rectMaxLat) {
        // 添加输入验证
        if (rectMinLon > rectMaxLon || rectMinLat > rectMaxLat) {
            throw new IllegalArgumentException("矩形参数无效：minLon不能大于maxLon，minLat不能大于maxLat");
        }

        return pointLon >= rectMinLon && pointLon <= rectMaxLon &&
                pointLat >= rectMinLat && pointLat <= rectMaxLat;
    }


    /**
     * 判断点是否在圆形内
     *
     * @param point   点坐标
     * @param center  圆心坐标
     * @param radiusM 圆的半径（米）
     *
     * @return 如果点在圆内则返回true，否则返回false
     */
    public boolean inCircle(CoordinatePoint point, CoordinatePoint center, double radiusM) {
        // 添加空值检查
        if (point == null || center == null) {
            throw new IllegalArgumentException("点坐标不能为null");
        }

        // 计算两点间距离
        double distance = haversine(point, center);
        return distance <= radiusM;
    }


    /**
     * 判断点是否在圆形内
     *
     * @param pointLon  点的经度
     * @param pointLat  点的纬度
     * @param centerLon 圆心经度
     * @param centerLat 圆心纬度
     * @param radiusM   圆的半径（米）
     *
     * @return 如果点在圆内则返回true，否则返回false
     */
    public boolean inCircle(double pointLon, double pointLat,
                            double centerLon, double centerLat,
                            double radiusM) {
        return inCircle(
                new CoordinatePoint(pointLon, pointLat),
                new CoordinatePoint(centerLon, centerLat),
                radiusM);
    }

    /**
     * 判断点是否在多边形内
     *
     * @param pointLon   点的经度
     * @param pointLat   点的纬度
     * @param polygonWkt 多边形的WKT表示
     *
     * @return 如果点在多边形内则返回true，否则返回false
     *
     * @throws Exception 如果WKT解析失败
     */
    public boolean inPolygon(double pointLon, double pointLat, String polygonWkt) throws Exception {
        if (polygonWkt == null || polygonWkt.isEmpty()) {
            throw new IllegalArgumentException("多边形WKT不能为null或空");
        }

        try {
            Geometry polygon = config.geometryFactory.createGeometry(new GeometryJSON().read(polygonWkt));
            if (polygon == null) {
                throw new IllegalArgumentException("解析多边形WKT失败");
            }

            Geometry point = config.geometryFactory.createPoint(new Coordinate(pointLon, pointLat));
            return point.within(polygon);
        } catch (Exception e) {
            throw new Exception("解析多边形WKT或判断点在多边形内时出错: " + e.getMessage(), e);
        }
    }


    /**
     * 创建缓冲区
     *
     * @param geom      原始几何形状
     * @param distance  缓冲距离（米）
     * @param originLon 原始经度，用于选择正确的投影坐标系
     *
     * @return 缓冲区几何形状
     *
     * @throws Exception 坐标转换异常
     */
    public Geometry buffer(Geometry geom, double distance, double originLon) throws Exception {
        try {
            CoordinateReferenceSystem src = CRS.decode(config.WGS84, true);
            CoordinateReferenceSystem tgt = pickCrs(originLon); // 使用传入的原始经度
            MathTransform tx = CRS.findMathTransform(src, tgt, true);
            Geometry proj = JTS.transform(geom, tx);
            Geometry buffer = proj.buffer(distance);
            Geometry wgs = JTS.transform(buffer, CRS.findMathTransform(tgt, src, true));
            return wgs;
        } catch (Exception e) {
            throw new Exception("创建缓冲区时出错: " + e.getMessage(), e);
        }
    }


    /**
     * 计算两个几何形状的交集
     *
     * @param g1 第一个几何形状
     * @param g2 第二个几何形状
     *
     * @return 两个几何形状的交集
     */
    public Geometry intersection(Geometry g1, Geometry g2) {
        return g1.intersection(g2);
    }

    /**
     * 计算两个几何形状的并集
     *
     * @param g1 第一个几何形状
     * @param g2 第二个几何形状
     *
     * @return 两个几何形状的并集
     */
    public Geometry union(Geometry g1, Geometry g2) {
        return g1.union(g2);
    }

    /**
     * 计算两个几何形状的差集
     *
     * @param g1 第一个几何形状
     * @param g2 第二个几何形状
     *
     * @return 两个几何形状的差集(g1-g2)
     */
    public Geometry difference(Geometry g1, Geometry g2) {
        return g1.difference(g2);
    }

    /**
     * 计算两个几何形状的对称差集
     *
     * @param g1 第一个几何形状
     * @param g2 第二个几何形状
     *
     * @return 两个几何形状的对称差集
     */
    public Geometry symDifference(Geometry g1, Geometry g2) {
        return g1.symDifference(g2);
    }

    /**
     * 将WGS84坐标转换为GCJ02坐标（火星坐标）
     *
     * @param lon WGS84经度
     * @param lat WGS84纬度
     *
     * @return GCJ02坐标点
     */
    public CoordinatePoint wgs84ToGcj02(double lon, double lat) {
        if (outOfChina(lon, lat)) {
            return new CoordinatePoint(lon, lat);
        }

        double[] delta = delta(lon, lat);
        return new CoordinatePoint(lon + delta[0], lat + delta[1]);
    }

    /**
     * 将GCJ02坐标转换为WGS84坐标
     *
     * @param lon GCJ02经度
     * @param lat GCJ02纬度
     *
     * @return WGS84坐标点
     */
    public CoordinatePoint gcj02ToWgs84(double lon, double lat) {
        if (outOfChina(lon, lat)) {
            return new CoordinatePoint(lon, lat);
        }

        double[] delta = delta(lon, lat);
        return new CoordinatePoint(lon - delta[0], lat - delta[1]);
    }

    /**
     * 将GCJ02坐标转换为BD09坐标（百度坐标）
     *
     * @param lon GCJ02经度
     * @param lat GCJ02纬度
     *
     * @return BD09坐标点
     */
    public CoordinatePoint gcj02ToBd09(double lon, double lat) {
        double z = Math.sqrt(lon * lon + lat * lat) + 0.00002 * Math.sin(lat * Math.PI);
        double theta = Math.atan2(lat, lon) + 0.000003 * Math.cos(lon * Math.PI);
        double bd_lon = z * Math.cos(theta) + 0.0065;
        double bd_lat = z * Math.sin(theta) + 0.006;
        return new CoordinatePoint(bd_lon, bd_lat);
    }

    /**
     * 将BD09坐标转换为GCJ02坐标
     *
     * @param lon BD09经度
     * @param lat BD09纬度
     *
     * @return GCJ02坐标点
     */
    public CoordinatePoint bd09ToGcj02(double lon, double lat) {
        double x = lon - 0.0065, y = lat - 0.006;
        double z = Math.sqrt(x * x + y * y) - 0.00002 * Math.sin(y * Math.PI);
        double theta = Math.atan2(y, x) - 0.000003 * Math.cos(x * Math.PI);
        double gg_lon = z * Math.cos(theta);
        double gg_lat = z * Math.sin(theta);
        return new CoordinatePoint(gg_lon, gg_lat);
    }

    /**
     * 将WGS84坐标直接转换为BD09坐标
     *
     * @param lon WGS84经度
     * @param lat WGS84纬度
     *
     * @return BD09坐标点
     */
    public CoordinatePoint wgs84ToBd09(double lon, double lat) {
        CoordinatePoint gcj02 = wgs84ToGcj02(lon, lat);
        return gcj02ToBd09(gcj02.getLon(), gcj02.getLat());
    }

    /**
     * 将BD09坐标直接转换为WGS84坐标
     *
     * @param lon BD09经度
     * @param lat BD09纬度
     *
     * @return WGS84坐标点
     */
    public CoordinatePoint bd09ToWgs84(double lon, double lat) {
        CoordinatePoint gcj02 = bd09ToGcj02(lon, lat);
        return gcj02ToWgs84(gcj02.getLon(), gcj02.getLat());
    }

    /**
     * 批量将WGS84坐标点列表转换为GCJ02坐标点列表
     *
     * @param points WGS84坐标点列表
     *
     * @return GCJ02坐标点列表
     */
    public List<CoordinatePoint> wgs84ToGcj02(List<CoordinatePoint> points) {
        return points.stream()
                .map(p -> wgs84ToGcj02(p.getLon(), p.getLat()))
                .collect(Collectors.toList());
    }

    /**
     * 批量将GCJ02坐标点列表转换为WGS84坐标点列表
     *
     * @param points GCJ02坐标点列表
     *
     * @return WGS84坐标点列表
     */
    public List<CoordinatePoint> gcj02ToWgs84(List<CoordinatePoint> points) {
        return points.stream()
                .map(p -> gcj02ToWgs84(p.getLon(), p.getLat()))
                .collect(Collectors.toList());
    }

    /**
     * 批量将GCJ02坐标点列表转换为BD09坐标点列表
     *
     * @param points GCJ02坐标点列表
     *
     * @return BD09坐标点列表
     */
    public List<CoordinatePoint> gcj02ToBd09(List<CoordinatePoint> points) {
        return points.stream()
                .map(p -> gcj02ToBd09(p.getLon(), p.getLat()))
                .collect(Collectors.toList());
    }

    /**
     * 批量将BD09坐标点列表转换为GCJ02坐标点列表
     *
     * @param points BD09坐标点列表
     *
     * @return GCJ02坐标点列表
     */
    public List<CoordinatePoint> bd09ToGcj02(List<CoordinatePoint> points) {
        return points.stream()
                .map(p -> bd09ToGcj02(p.getLon(), p.getLat()))
                .collect(Collectors.toList());
    }

    /**
     * 批量将WGS84坐标点列表直接转换为BD09坐标点列表
     *
     * @param points WGS84坐标点列表
     *
     * @return BD09坐标点列表
     */
    public List<CoordinatePoint> wgs84ToBd09(List<CoordinatePoint> points) {
        return points.stream()
                .map(p -> wgs84ToBd09(p.getLon(), p.getLat()))
                .collect(Collectors.toList());
    }

    /**
     * 批量将BD09坐标点列表直接转换为WGS84坐标点列表
     *
     * @param points BD09坐标点列表
     *
     * @return WGS84坐标点列表
     */
    public List<CoordinatePoint> bd09ToWgs84(List<CoordinatePoint> points) {
        return points.stream()
                .map(p -> bd09ToWgs84(p.getLon(), p.getLat()))
                .collect(Collectors.toList());
    }

    /**
     * 批量将TrackPoint列表从WGS84坐标转换为GCJ02坐标（修改原对象）
     *
     * @param points WGS84坐标的TrackPoint列表
     *
     * @return WGS84坐标的TrackPoint列表（原列表，坐标已更新）
     */
    public List<TrackPoint> wgs84ToGcj02TrackPoints(List<TrackPoint> points) {
        return points.stream()
                .peek(p -> {
                    CoordinatePoint result = wgs84ToGcj02(p.getLon(), p.getLat());
                    // 更新现有TrackPoint对象的坐标
                    p.setLon(result.getLon());
                    p.setLat(result.getLat());
                })
                .collect(Collectors.toList());
    }

    /**
     * 批量将TrackPoint列表从GCJ02坐标转换为WGS84坐标（修改原对象）
     *
     * @param points GCJ02坐标的TrackPoint列表
     *
     * @return WGS84坐标的TrackPoint列表（原列表，坐标已更新）
     */
    public List<TrackPoint> gcj02ToWgs84TrackPoints(List<TrackPoint> points) {
        return points.stream()
                .peek(p -> {
                    CoordinatePoint result = gcj02ToWgs84(p.getLon(), p.getLat());
                    // 更新现有TrackPoint对象的坐标
                    p.setLon(result.getLon());
                    p.setLat(result.getLat());
                })
                .collect(Collectors.toList());
    }

    /**
     * 批量将TrackPoint列表从GCJ02坐标转换为BD09坐标（创建新对象）
     *
     * @param points GCJ02坐标的TrackPoint列表
     *
     * @return BD09坐标的TrackPoint列表（新列表）
     */
    public List<TrackPoint> gcj02ToBd09TrackPoints(List<TrackPoint> points) {
        return points.stream()
                .map(p -> {
                    CoordinatePoint result = gcj02ToBd09(p.getLon(), p.getLat());
                    TrackPoint tp = new TrackPoint(result.getLon(), result.getLat());
                    BeanUtil.copyProperties(p, tp, "lon", "lat");
                    return tp;
                })
                .collect(Collectors.toList());
    }

    /**
     * 批量将TrackPoint列表从BD09坐标转换为GCJ02坐标（创建新对象）
     *
     * @param points BD09坐标的TrackPoint列表
     *
     * @return GCJ02坐标的TrackPoint列表（新列表）
     */
    public List<TrackPoint> bd09ToGcj02TrackPoints(List<TrackPoint> points) {
        return points.stream()
                .map(p -> {
                    CoordinatePoint result = gcj02ToBd09(p.getLon(), p.getLat());
                    TrackPoint tp = new TrackPoint(result.getLon(), result.getLat());
                    BeanUtil.copyProperties(p, tp, "lon", "lat");
                    return tp;
                })
                .collect(Collectors.toList());
    }

    /**
     * 批量将TrackPoint列表从WGS84坐标直接转换为BD09坐标（创建新对象）
     *
     * @param points WGS84坐标的TrackPoint列表
     *
     * @return BD09坐标的TrackPoint列表（新列表）
     */
    public List<TrackPoint> wgs84ToBd09TrackPoints(List<TrackPoint> points) {
        return points.stream()
                .map(p -> {
                    CoordinatePoint result = gcj02ToBd09(p.getLon(), p.getLat());
                    TrackPoint tp = new TrackPoint(result.getLon(), result.getLat());
                    BeanUtil.copyProperties(p, tp, "lon", "lat");
                    return tp;
                })
                .collect(Collectors.toList());
    }

    /**
     * 批量将TrackPoint列表从BD09坐标直接转换为WGS84坐标（创建新对象）
     *
     * @param points BD09坐标的TrackPoint列表
     *
     * @return WGS84坐标的TrackPoint列表（新列表）
     */
    public List<TrackPoint> bd09ToWgs84TrackPoints(List<TrackPoint> points) {
        return points.stream()
                .map(p -> {
                    CoordinatePoint result = gcj02ToBd09(p.getLon(), p.getLat());
                    TrackPoint tp = new TrackPoint(result.getLon(), result.getLat());
                    BeanUtil.copyProperties(p, tp, "lon", "lat");
                    return tp;
                })
                .collect(Collectors.toList());
    }


    /**
     * 判断点是否在中国境外
     * <p>
     * 中国范围定义为:
     * 经度: 72.004 < lon < 137.8347
     * 纬度: 0.8293 < lat < 55.8271
     *
     * @param lon 经度
     * @param lat 纬度
     *
     * @return 如果在中国境外返回true，否则返回false
     */
    public boolean outOfChina(double lon, double lat) {
        // 严格边界检查，边界上的点被认为在中国境内
        return lon < 72.004 || lon > 137.8347 || lat < 0.8293 || lat > 55.8271;
    }


}
