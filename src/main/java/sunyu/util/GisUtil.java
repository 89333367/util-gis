package sunyu.util;

import cn.hutool.core.bean.BeanUtil;
import cn.hutool.log.Log;
import cn.hutool.log.LogFactory;
import org.geotools.geojson.geom.GeometryJSON;
import org.geotools.geometry.jts.JTS;
import org.geotools.referencing.CRS;
import org.locationtech.jts.geom.*;
import org.locationtech.jts.io.WKTReader;
import org.opengis.referencing.crs.CoordinateReferenceSystem;
import org.opengis.referencing.operation.MathTransform;
import sunyu.util.pojo.CoordinatePoint;
import sunyu.util.pojo.TrackPoint;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.concurrent.ConcurrentHashMap;
import java.util.stream.Collectors;

/**
 * GIS工具类，用于轨迹处理、空间计算等
 * 使用统一的Web Mercator投影系统（EPSG:3857）简化坐标转换
 *
 * @author SunYu
 */
public class GisUtil implements AutoCloseable {
    // 日志记录器，用于记录工具类的运行状态和调试信息
    private final Log log = LogFactory.get();
    // 配置参数，包含各种常量和默认值
    private final Config config;
    // 几何图形合并工具
    private final GeometryUnionUtil geometryUnionUtil;

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
        this.geometryUnionUtil = new GeometryUnionUtil(config.geometryFactory);
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

        // Web Mercator投影系统的EPSG代码，用于统一投影坐标系统
        private final String WEB_MERCATOR = "EPSG:3857";

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
     * @param code EPSG代码
     *
     * @return CoordinateReferenceSystem对象
     */
    private CoordinateReferenceSystem getCachedCRS(String code) {
        return config.crsCache.computeIfAbsent(code, key -> {
            try {
                return CRS.decode(key, true);
            } catch (Exception e) {
                throw new RuntimeException("Failed to parse CRS: " + key, e);
            }
        });
    }

    /**
     * 将WGS84坐标转换为Web Mercator投影坐标
     *
     * @param g WGS84坐标系下的几何图形
     *
     * @return Web Mercator投影坐标系下的几何图形
     *
     * @throws Exception 坐标转换异常
     */
    private Geometry wgs84ToWebMercator(Geometry g) throws Exception {
        CoordinateReferenceSystem src = getCachedCRS(config.WGS84);
        CoordinateReferenceSystem tgt = getCachedCRS(config.WEB_MERCATOR);
        MathTransform tx = CRS.findMathTransform(src, tgt, true);
        Geometry result = JTS.transform(g, tx);
        return result;
    }

    /**
     * 将Web Mercator投影坐标转换为WGS84坐标
     *
     * @param g Web Mercator投影坐标系下的几何图形
     *
     * @return WGS84坐标系下的几何图形
     *
     * @throws Exception 坐标转换异常
     */
    private Geometry webMercatorToWgs84(Geometry g) throws Exception {
        try {
            CoordinateReferenceSystem src = getCachedCRS(config.WEB_MERCATOR);
            CoordinateReferenceSystem tgt = getCachedCRS(config.WGS84);
            MathTransform tx = CRS.findMathTransform(src, tgt, true);
            Geometry result = JTS.transform(g, tx);

            // 验证转换后的坐标是否在合理范围内
            if (!isValidWgs84Geometry(result)) {
                // 如果转换失败，尝试使用不同的转换方法
                MathTransform tx2 = CRS.findMathTransform(src, tgt, false);
                result = JTS.transform(g, tx2);
            }

            return result;
        } catch (Exception e) {
            log.error("Web Mercator到WGS84坐标转换失败: " + e.getMessage(), e);
            throw e;
        }
    }


    /**
     * 验证WGS84坐标系下的几何图形坐标是否在合理范围内
     *
     * @param g WGS84坐标系下的几何图形
     *
     * @return 如果坐标在合理范围内返回true，否则返回false
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
     * 构建轨迹轮廓（重载方法，使用统一宽度）
     * 通过轨迹点生成带宽度的轮廓多边形，左右宽度相同（各为总宽度的一半）
     *
     * @param seg         轨迹点列表，至少需要3个点，要求点位坐标系为WGS84
     * @param totalWidthM 总宽度（米），轨迹线两侧的扩展距离之和
     *
     * @return Geometry对象，表示生成的轮廓多边形
     *
     * @throws Exception 可能抛出异常，如坐标转换错误或几何操作失败
     */
    public Geometry buildOutline(List<TrackPoint> seg, double totalWidthM) throws Exception {
        double halfWidth = totalWidthM / 2.0;
        return buildOutline(seg, halfWidth, halfWidth, 5.0, 20.0);
    }

    /**
     * 构建轨迹轮廓（重载方法，使用统一宽度）
     * 通过轨迹点生成带宽度的轮廓多边形，左右宽度相同（各为总宽度的一半）
     *
     * @param seg            轨迹点列表，至少需要3个点，要求点位坐标系为WGS84
     * @param totalWidthM    总宽度（米），轨迹线两侧的扩展距离之和
     * @param simplifyM      简化阈值（米），用于轨迹点简化以提高性能
     * @param maxEdgeLengthM 最大边缘长度（米），控制轮廓的精细程度
     *
     * @return Geometry对象，表示生成的轮廓多边形
     *
     * @throws Exception 可能抛出异常，如坐标转换错误或几何操作失败
     */
    public Geometry buildOutline(List<TrackPoint> seg,
                                 double totalWidthM,
                                 double simplifyM,
                                 double maxEdgeLengthM) throws Exception {
        double halfWidth = totalWidthM / 2.0;
        return buildOutline(seg, halfWidth, halfWidth, simplifyM, maxEdgeLengthM);
    }

    /**
     * 构建轨迹轮廓
     * 通过轨迹点生成带宽度的轮廓多边形，用于面积计算等操作
     * 使用统一的Web Mercator投影系统进行计算，避免复杂的投影选择
     *
     * @param seg            轨迹点列表，至少需要3个点，坐标系必须为WGS84
     * @param leftWidthM     左侧宽度（米），轨迹线左侧的扩展距离
     * @param rightWidthM    右侧宽度（米），轨迹线右侧的扩展距离
     * @param simplifyM      简化阈值（米），用于轨迹点简化以提高性能
     * @param maxEdgeLengthM 最大边缘长度（米），控制轮廓的精细程度
     *
     * @return Geometry对象，表示生成的轮廓多边形，返回WGS84坐标系中的几何对象
     *
     * @throws Exception 可能抛出异常，如坐标转换错误或几何操作失败
     */
    public Geometry buildOutline(List<TrackPoint> seg,
                                 double leftWidthM,
                                 double rightWidthM,
                                 double simplifyM,
                                 double maxEdgeLengthM) throws Exception {
        long startTime = System.currentTimeMillis();
        log.debug("[buildOutline] 开始处理轨迹轮廓，轨迹点数: {}, 左侧宽度: {} 米, 右侧宽度: {} 米, 简化阈值: {} 米, 最大边缘长度: {} 米",
                seg != null ? seg.size() : 0, leftWidthM, rightWidthM, simplifyM, maxEdgeLengthM);

        // 如果轨迹点少于3个，无法构成有效几何形状，抛出异常
        if (seg == null) {
            throw new IllegalArgumentException("轨迹段至少需要3个点");
        }
        if (leftWidthM < 0 || rightWidthM < 0 || simplifyM < 0 || maxEdgeLengthM < 0) {
            throw new IllegalArgumentException("所有参数必须为非负数");
        }

        // 使用Kalman滤波对轨迹点进行预处理，消除噪声影响
        long filterStartTime = System.currentTimeMillis();
        KalmanFilterUtil kalmanFilter = new KalmanFilterUtil();
        List<TrackPoint> filteredSeg = kalmanFilter.filterTrack(seg);
        long filterTime = System.currentTimeMillis() - filterStartTime;
        log.debug("[buildOutline] Kalman滤波处理完成，原始点数: {}, 滤波后点数: {}, 滤波耗时: {}ms",
                seg.size(), filteredSeg.size(), filterTime);

        // 先按时间戳排序，确保轨迹点的时序正确
        // 过滤掉无效的轨迹点（经纬度为0或超出范围的点）
        long processStartTime = System.currentTimeMillis();
        List<TrackPoint> sortedSeg = filteredSeg.stream()
                .sorted(Comparator.comparing(TrackPoint::getTime))
                .filter(p -> Math.abs(p.getLon()) <= 180 && Math.abs(p.getLat()) <= 90)
                .filter(p -> !(p.getLon() == 0 && p.getLat() == 0)) // 过滤掉经纬度都为0的点
                .collect(Collectors.toList());
        long processTime = System.currentTimeMillis() - processStartTime;

        log.debug("[buildOutline] 轨迹点过滤完成，滤波后点数: {}, 过滤后点数: {}, 过滤耗时: {}ms",
                filteredSeg.size(), sortedSeg.size(), processTime);

        if (sortedSeg.size() < 3) {
            throw new IllegalArgumentException("轨迹段至少需要3个有效点");
        }

        try {
            double minLon = sortedSeg.stream().mapToDouble(TrackPoint::getLon).min().orElse(0);
            double maxLon = sortedSeg.stream().mapToDouble(TrackPoint::getLon).max().orElse(0);
            double minLat = sortedSeg.stream().mapToDouble(TrackPoint::getLat).min().orElse(0);
            double maxLat = sortedSeg.stream().mapToDouble(TrackPoint::getLat).max().orElse(0);

            log.debug("[buildOutline] 轨迹范围: 经度[{}, {}], 纬度[{}, {}]", minLon, maxLon, minLat, maxLat);

            // 使用新的基于点缓冲区的方法
            long buildStartTime = System.currentTimeMillis();
            Geometry result = buildOutlineByPointBuffers(sortedSeg, leftWidthM, rightWidthM);
            long buildTime = System.currentTimeMillis() - buildStartTime;

            log.debug("[buildOutline] 轮廓构建完成，构建耗时: {}ms", buildTime);

            // 如果结果是MultiPolygon且包含多个部分，直接返回
            if (result instanceof MultiPolygon && result.getNumGeometries() > 1) {
                log.debug("[buildOutline] 结果为MultiPolygon，包含 {} 个部分", result.getNumGeometries());
                long endTime = System.currentTimeMillis();
                log.debug("[buildOutline] 总计耗时: {}ms", endTime - startTime);
                return result;
            }

            // 如果是单个Polygon，检查是否需要进一步分割
            if (result instanceof Polygon || (result instanceof MultiPolygon && result.getNumGeometries() == 1)) {
                log.debug("[buildOutline] 结果为Polygon或单个MultiPolygon");
                long endTime = System.currentTimeMillis();
                log.debug("[buildOutline] 总计耗时: {}ms", endTime - startTime);
                // 可以在这里添加额外的处理逻辑
                return result;
            }

            long endTime = System.currentTimeMillis();
            log.debug("[buildOutline] 处理完成，总计耗时: {}ms", endTime - startTime);

            return result;
        } catch (Exception e) {
            log.error("构建轨迹轮廓失败: " + e.getMessage(), e);
            throw new Exception("构建轨迹轮廓失败: " + e.getMessage(), e);
        }
    }

    /**
     * 基于点缓冲区构建轨迹轮廓
     * 为每个轨迹点创建缓冲区，然后合并这些缓冲区形成最终的轮廓
     *
     * @param seg         轨迹点列表
     * @param leftWidthM  左侧宽度（米）
     * @param rightWidthM 右侧宽度（米）
     *
     * @return 生成的几何对象
     *
     * @throws Exception 坐标转换异常
     */
    private Geometry buildOutlineByPointBuffers(List<TrackPoint> seg, double leftWidthM, double rightWidthM) throws Exception {
        long startTime = System.currentTimeMillis();
        try {
            // 如果左右宽度相同，使用简单的点缓冲区方法
            if (Math.abs(leftWidthM - rightWidthM) < 0.001) {
                log.debug("[buildOutlineByPointBuffers] 使用简单缓冲区方法 (左右宽度相同)");
                Geometry result = buildOutlineBySimpleBuffers(seg, leftWidthM);
                long endTime = System.currentTimeMillis();
                log.debug("[buildOutlineByPointBuffers] 简单缓冲区方法完成，耗时: {}ms", endTime - startTime);
                return result;
            } else {
                // 如果左右宽度不同，使用更复杂的方法
                log.debug("[buildOutlineByPointBuffers] 使用复杂缓冲区方法 (左右宽度不同)");
                Geometry result = buildOutlineByComplexBuffers(seg, leftWidthM, rightWidthM);
                long endTime = System.currentTimeMillis();
                log.debug("[buildOutlineByPointBuffers] 复杂缓冲区方法完成，耗时: {}ms", endTime - startTime);
                return result;
            }
        } catch (Exception e) {
            log.error("基于点缓冲区构建轮廓失败: " + e.getMessage(), e);
            throw e;
        }
    }

    /**
     * 使用相同左右宽度构建轮廓
     *
     * @param seg    轨迹点列表
     * @param widthM 宽度（米）
     *
     * @return 生成的几何对象
     *
     * @throws Exception 坐标转换异常
     */
    private Geometry buildOutlineBySimpleBuffers(List<TrackPoint> seg, double widthM) throws Exception {
        long startTime = System.currentTimeMillis();
        log.debug("[buildOutlineBySimpleBuffers] 开始处理 {} 个轨迹点，缓冲区宽度: {} 米", seg.size(), widthM);

        // 为每个点创建缓冲区
        List<Geometry> pointBuffers = new ArrayList<>();

        long convertTime = 0;
        long bufferTime = 0;
        long totalTime = 0;

        for (int i = 0; i < seg.size(); i++) {
            long startPoint = System.currentTimeMillis();

            TrackPoint point = seg.get(i);
            Coordinate coord = new Coordinate(point.getLon(), point.getLat());
            Geometry pointGeom = config.geometryFactory.createPoint(coord);

            // 转换到Web Mercator投影坐标系
            long startConvert = System.currentTimeMillis();
            Geometry projPoint = wgs84ToWebMercator(pointGeom);
            convertTime += System.currentTimeMillis() - startConvert;

            // 创建缓冲区
            long startBuffer = System.currentTimeMillis();
            Geometry buffer = projPoint.buffer(widthM);
            bufferTime += System.currentTimeMillis() - startBuffer;

            pointBuffers.add(buffer);

            totalTime += System.currentTimeMillis() - startPoint;

            // 每处理100个点打印一次进度
            if (i > 0 && i % 100 == 0) {
                log.trace("[buildOutlineBySimpleBuffers] 已处理 {}/{} 个点, 转换耗时: {}ms, 缓冲耗时: {}ms, 平均每个点耗时: {}ms",
                        i, seg.size(), convertTime, bufferTime, totalTime / i);
            }
        }

        log.debug("[buildOutlineBySimpleBuffers] 点缓冲区创建完成，共 {} 个缓冲区, 转换总耗时: {}ms, 缓冲总耗时: {}ms",
                pointBuffers.size(), convertTime, bufferTime);

        // 使用GeometryCollection优化合并所有缓冲区
        long startUnion = System.currentTimeMillis();
        log.debug("[buildOutlineBySimpleBuffers] 开始使用GeometryCollection优化合并 {} 个缓冲区", pointBuffers.size());

        // 使用优化的合并方法
        Geometry union = geometryUnionUtil.union(pointBuffers);

        long unionTime = System.currentTimeMillis() - startUnion;
        log.debug("[buildOutlineBySimpleBuffers] 缓冲区合并完成，合并耗时: {}ms", unionTime);

        // 转换回WGS84坐标系
        long startBackConvert = System.currentTimeMillis();
        Geometry result = webMercatorToWgs84(union);
        long backConvertTime = System.currentTimeMillis() - startBackConvert;

        long endTime = System.currentTimeMillis();
        log.debug("[buildOutlineBySimpleBuffers] 处理完成，总计耗时: {}ms (转换:{}ms, 缓冲:{}ms, 合并:{}ms, 回转:{}ms)",
                endTime - startTime, convertTime, bufferTime, unionTime, backConvertTime);

        return result;
    }

    /**
     * 使用不同左右宽度构建轮廓
     *
     * @param seg         轨迹点列表
     * @param leftWidthM  左侧宽度（米）
     * @param rightWidthM 右侧宽度（米）
     *
     * @return 生成的几何对象
     *
     * @throws Exception 坐标转换异常
     */
    private Geometry buildOutlineByComplexBuffers(List<TrackPoint> seg, double leftWidthM, double rightWidthM) throws Exception {
        long startTime = System.currentTimeMillis();
        log.debug("[buildOutlineByComplexBuffers] 开始处理 {} 个轨迹点，左侧宽度: {} 米，右侧宽度: {} 米", seg.size(), leftWidthM, rightWidthM);

        // 为每个线段创建缓冲区
        List<Geometry> segmentBuffers = new ArrayList<>();

        long convertTime = 0;
        long bufferTime = 0;
        long totalTime = 0;

        for (int i = 0; i < seg.size() - 1; i++) {
            long startPoint = System.currentTimeMillis();

            TrackPoint point1 = seg.get(i);
            TrackPoint point2 = seg.get(i + 1);

            // 创建线段
            Coordinate[] coords = {
                    new Coordinate(point1.getLon(), point1.getLat()),
                    new Coordinate(point2.getLon(), point2.getLat())
            };
            LineString line = config.geometryFactory.createLineString(coords);

            // 转换到Web Mercator投影坐标系
            long startConvert = System.currentTimeMillis();
            LineString projLine = (LineString) wgs84ToWebMercator(line);
            convertTime += System.currentTimeMillis() - startConvert;

            // 创建缓冲区（使用左右宽度的平均值作为近似）
            long startBuffer = System.currentTimeMillis();
            double avgWidth = (leftWidthM + rightWidthM) / 2.0;
            Geometry buffer = projLine.buffer(avgWidth);
            bufferTime += System.currentTimeMillis() - startBuffer;

            segmentBuffers.add(buffer);

            totalTime += System.currentTimeMillis() - startPoint;

            // 每处理100个点打印一次进度
            if (i > 0 && i % 100 == 0) {
                log.trace("[buildOutlineByComplexBuffers] 已处理 {}/{} 个线段, 转换耗时: {}ms, 缓冲耗时: {}ms, 平均每个线段耗时: {}ms",
                        i, seg.size() - 1, convertTime, bufferTime, totalTime / i);
            }
        }

        log.debug("[buildOutlineByComplexBuffers] 线段缓冲区创建完成，共 {} 个缓冲区, 转换总耗时: {}ms, 缓冲总耗时: {}ms",
                segmentBuffers.size(), convertTime, bufferTime);

        // 使用GeometryCollection优化合并所有缓冲区
        long startUnion = System.currentTimeMillis();
        log.debug("[buildOutlineByComplexBuffers] 开始使用GeometryCollection优化合并 {} 个缓冲区", segmentBuffers.size());

        // 使用优化的合并方法
        Geometry union = geometryUnionUtil.union(segmentBuffers);

        long unionTime = System.currentTimeMillis() - startUnion;
        log.debug("[buildOutlineByComplexBuffers] 缓冲区合并完成，合并耗时: {}ms", unionTime);

        // 转换回WGS84坐标系
        long startBackConvert = System.currentTimeMillis();
        Geometry result = webMercatorToWgs84(union);
        long backConvertTime = System.currentTimeMillis() - startBackConvert;

        long endTime = System.currentTimeMillis();
        log.debug("[buildOutlineByComplexBuffers] 处理完成，总计耗时: {}ms (转换:{}ms, 缓冲:{}ms, 合并:{}ms, 回转:{}ms)",
                endTime - startTime, convertTime, bufferTime, unionTime, backConvertTime);

        return result;
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
        double dLat = Math.toRadians(p2.getLat() - p1.getLat());
        double dLon = Math.toRadians(p2.getLon() - p1.getLon());
        double a = Math.sin(dLat / 2) * Math.sin(dLat / 2) +
                Math.cos(Math.toRadians(p1.getLat())) *
                        Math.cos(Math.toRadians(p2.getLat())) *
                        Math.sin(dLon / 2) * Math.sin(dLon / 2);
        return config.R * 2 * Math.atan2(Math.sqrt(a), Math.sqrt(1 - a));
    }


    /**
     * 合并相交的多边形
     *
     * @param polygons 多边形列表
     *
     * @return 合并后的多边形列表
     */
    private List<Geometry> mergeIntersectingPolygons(List<Geometry> polygons) {
        List<Geometry> result = new ArrayList<>();

        for (Geometry polygon : polygons) {
            boolean merged = false;

            for (int i = 0; i < result.size(); i++) {
                Geometry existingPolygon = result.get(i);

                // 检查是否相交或接触
                if (existingPolygon.intersects(polygon) || existingPolygon.touches(polygon)) {
                    // 合并多边形
                    Geometry union = existingPolygon.union(polygon);
                    result.set(i, union);
                    merged = true;
                    break;
                }
            }

            if (!merged) {
                result.add(polygon);
            }
        }

        // 递归合并直到没有更多可以合并的多边形
        if (result.size() < polygons.size()) {
            return mergeIntersectingPolygons(result);
        }

        return result;
    }


    /**
     * 将几何图形转换为WKT格式
     * 确保输出为WGS84坐标系的标准地理坐标
     *
     * @param g 几何图形对象
     *
     * @return WKT字符串，表示几何图形
     */
    public String toWkt(Geometry g) {
        try {
            // 检查几何图形是否已经是WGS84坐标系并且坐标有效
            if (isLikelyWgs84(g) && hasValidCoordinates(g)) {
                String wkt = g.toText();
                return wkt;
            }

            // 转换到WGS84坐标系
            Geometry wgs = webMercatorToWgs84(g);
            // 验证转换后的几何图形
            if (!isValidWgs84Geometry(wgs) || !hasValidCoordinates(wgs)) {
                // 尝试清理几何图形
                wgs = wgs.buffer(0);
                if (!isValidWgs84Geometry(wgs) || !hasValidCoordinates(wgs)) {
                    String wkt = g.toText();
                    return wkt;
                }
            }

            String wkt = wgs.toText();
            return wkt;
        } catch (Exception e) {
            log.error("坐标转换失败: " + e.getMessage(), e);
            String wkt = g.toText();
            return wkt;
        }
    }

    /**
     * 检查几何图形中的坐标是否有效（没有接近零的异常值）
     *
     * @param g 几何图形对象
     *
     * @return 如果坐标有效返回true，否则返回false
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
     * 判断几何图形是否可能已经是WGS84坐标系
     *
     * @param g 几何图形对象
     *
     * @return 如果可能是WGS84坐标系返回true，否则返回false
     */
    private boolean isLikelyWgs84(Geometry g) {
        Coordinate[] coordinates = g.getCoordinates();
        int invalidCount = 0;
        int totalCount = coordinates.length;

        for (Coordinate coord : coordinates) {
            // 更严格的判断：如果坐标值非常接近0（在误差范围内）或者超出地理范围，则认为不是WGS84坐标系
            if (Math.abs(coord.x) > 180 || Math.abs(coord.y) > 90 ||
                    (Math.abs(coord.x) < 1e-10 && Math.abs(coord.x) > 0) ||
                    (Math.abs(coord.y) < 1e-10 && Math.abs(coord.y) > 0) ||
                    (coord.x == 0 && coord.y == 0)) {
                invalidCount++;
            }
        }

        // 如果超过5%的坐标点有问题，认为不是WGS84坐标系
        return (double) invalidCount / totalCount < 0.05;
    }


    /**
     * 计算几何图形的面积（mu单位）
     * 通过坐标转换将几何图形转换为投影坐标后计算面积，结果以亩为单位
     *
     * @param outline 几何图形
     *
     * @return 面积（mu），以亩为单位，保留3位小数
     *
     * @throws RuntimeException 如果坐标转换过程中发生错误
     */
    public double calcMu(Geometry outline) throws RuntimeException {
        try {
            // 转换到Web Mercator投影坐标系计算面积
            Geometry proj = wgs84ToWebMercator(outline);

            // 计算面积并转换为mu单位（亩）
            return Math.round(proj.getArea() * config.MU_PER_SQ_METER * 1000.0) / 1000.0;
        } catch (Exception e) {
            throw new RuntimeException("计算面积时出错: " + e.getMessage(), e);
        }
    }

    /**
     * wkt计算亩数（WGS84坐标系）
     * 直接计算WGS84坐标系下几何图形的面积，结果以亩为单位
     *
     * @param wkt wkt字符串，要求为WGS84坐标系
     *
     * @return 面积（mu），以亩为单位，保留3位小数
     */
    public double calcMu(String wkt) {
        Geometry geometry = fromWkt(wkt);
        return calcMu(geometry);
    }

    /**
     * 从WKT字符串创建Geometry对象,WKT字符串必须是WGS84坐标系
     * 支持POLYGON和MULTIPOLYGON两种几何类型
     *
     * @param wkt WKT字符串，必须是POLYGON或MULTIPOLYGON类型
     *
     * @return 解析后的Geometry对象
     *
     * @throws IllegalArgumentException 如果WKT字符串为空或不支持的几何类型
     * @throws RuntimeException         如果解析过程中发生错误
     */
    public Geometry fromWkt(String wkt) {
        if (wkt == null || wkt.trim().isEmpty()) {
            throw new IllegalArgumentException("WKT字符串不能为空");
        }

        try {
            WKTReader wktReader = new WKTReader(config.geometryFactory);
            Geometry geometry = wktReader.read(wkt);

            if (geometry == null) {
                throw new RuntimeException("无法解析WKT字符串: " + wkt);
            }

            String geometryType = geometry.getGeometryType().toUpperCase();
            if (!("POLYGON".equals(geometryType) || "MULTIPOLYGON".equals(geometryType))) {
                throw new IllegalArgumentException("不支持的几何类型: " + geometryType + "。仅支持POLYGON和MULTIPOLYGON");
            }

            return geometry;
        } catch (org.locationtech.jts.io.ParseException e) {
            throw new RuntimeException("解析WKT字符串时出错: " + e.getMessage(), e);
        }
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
        double latRad = Math.toRadians(lat);
        double metersPerDegree = 2 * Math.PI * config.R * Math.cos(latRad) / 360;

        if (Math.abs(metersPerDegree) < 1e-10) {
            return 1.0 / (2 * Math.PI * config.R / 360);
        }

        return 1.0 / metersPerDegree;
    }

    /**
     * 判断两个几何形状是否拓扑相等
     */
    public boolean equals(Geometry g1, Geometry g2) {
        return g1.equals(g2);
    }

    /**
     * 判断两个几何形状是否没有共有点（脱节）
     */
    public boolean disjoint(Geometry g1, Geometry g2) {
        return g1.disjoint(g2);
    }

    /**
     * 判断两个几何形状是否至少有一个共有点（相交）
     */
    public boolean intersects(Geometry g1, Geometry g2) {
        return g1.intersects(g2);
    }

    /**
     * 判断两个几何形状是否有至少一个公共的边界点，但是没有内部点（接触）
     */
    public boolean touches(Geometry g1, Geometry g2) {
        return g1.touches(g2);
    }

    /**
     * 判断两个几何形状是否共享一些但不是所有的内部点（交叉）
     */
    public boolean crosses(Geometry g1, Geometry g2) {
        return g1.crosses(g2);
    }

    /**
     * 判断几何形状A是否完全在几何形状B内部（内含）
     */
    public boolean within(Geometry g1, Geometry g2) {
        return g1.within(g2);
    }

    /**
     * 判断几何形状A是否包含几何形状B（包含）
     */
    public boolean contains(Geometry g1, Geometry g2) {
        return g1.contains(g2);
    }

    /**
     * 判断两个几何形状是否共享一部分但不是所有的公共点，而且相交处有他们自己相同的区域（重叠）
     */
    public boolean overlaps(Geometry g1, Geometry g2) {
        return g1.overlaps(g2);
    }

    /**
     * 判断点是否在矩形内
     */
    public boolean inRectangle(double pointLon, double pointLat,
                               double rectMinLon, double rectMinLat,
                               double rectMaxLon, double rectMaxLat) {
        if (rectMinLon > rectMaxLon || rectMinLat > rectMaxLat) {
            throw new IllegalArgumentException("矩形参数无效：minLon不能大于maxLon，minLat不能大于maxLat");
        }

        return pointLon >= rectMinLon && pointLon <= rectMaxLon &&
                pointLat >= rectMinLat && pointLat <= rectMaxLat;
    }

    /**
     * 判断点是否在圆形内
     */
    public boolean inCircle(CoordinatePoint point, CoordinatePoint center, double radiusM) {
        if (point == null || center == null) {
            throw new IllegalArgumentException("点坐标不能为null");
        }

        double distance = haversine(point, center);
        return distance <= radiusM;
    }

    /**
     * 判断点是否在圆形内
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
     */
    public Geometry buffer(Geometry geom, double distance, double originLon) throws Exception {
        try {
            // 转换到Web Mercator投影坐标系
            Geometry proj = wgs84ToWebMercator(geom);
            Geometry buffer = proj.buffer(distance);
            // 转换回WGS84坐标系
            return webMercatorToWgs84(buffer);
        } catch (Exception e) {
            throw new Exception("创建缓冲区时出错: " + e.getMessage(), e);
        }
    }

    /**
     * 计算两个几何形状的交集
     */
    public Geometry intersection(Geometry g1, Geometry g2) {
        return g1.intersection(g2);
    }

    /**
     * 计算两个几何形状的并集
     */
    public Geometry union(Geometry g1, Geometry g2) {
        return g1.union(g2);
    }

    /**
     * 计算两个几何形状的差集
     */
    public Geometry difference(Geometry g1, Geometry g2) {
        return g1.difference(g2);
    }

    /**
     * 计算两个几何形状的对称差集
     */
    public Geometry symDifference(Geometry g1, Geometry g2) {
        return g1.symDifference(g2);
    }

    /**
     * 将WGS84坐标转换为GCJ02坐标（火星坐标）
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
     */
    public CoordinatePoint wgs84ToBd09(double lon, double lat) {
        CoordinatePoint gcj02 = wgs84ToGcj02(lon, lat);
        return gcj02ToBd09(gcj02.getLon(), gcj02.getLat());
    }

    /**
     * 将BD09坐标直接转换为WGS84坐标
     */
    public CoordinatePoint bd09ToWgs84(double lon, double lat) {
        CoordinatePoint gcj02 = bd09ToGcj02(lon, lat);
        return gcj02ToWgs84(gcj02.getLon(), gcj02.getLat());
    }

    /**
     * 批量将WGS84坐标点列表转换为GCJ02坐标点列表
     */
    public List<CoordinatePoint> wgs84ToGcj02(List<CoordinatePoint> points) {
        return points.stream()
                .map(p -> wgs84ToGcj02(p.getLon(), p.getLat()))
                .collect(Collectors.toList());
    }

    /**
     * 批量将GCJ02坐标点列表转换为WGS84坐标点列表
     */
    public List<CoordinatePoint> gcj02ToWgs84(List<CoordinatePoint> points) {
        return points.stream()
                .map(p -> gcj02ToWgs84(p.getLon(), p.getLat()))
                .collect(Collectors.toList());
    }

    /**
     * 批量将GCJ02坐标点列表转换为BD09坐标点列表
     */
    public List<CoordinatePoint> gcj02ToBd09(List<CoordinatePoint> points) {
        return points.stream()
                .map(p -> gcj02ToBd09(p.getLon(), p.getLat()))
                .collect(Collectors.toList());
    }

    /**
     * 批量将BD09坐标点列表转换为GCJ02坐标点列表
     */
    public List<CoordinatePoint> bd09ToGcj02(List<CoordinatePoint> points) {
        return points.stream()
                .map(p -> bd09ToGcj02(p.getLon(), p.getLat()))
                .collect(Collectors.toList());
    }

    /**
     * 批量将WGS84坐标点列表直接转换为BD09坐标点列表
     */
    public List<CoordinatePoint> wgs84ToBd09(List<CoordinatePoint> points) {
        return points.stream()
                .map(p -> wgs84ToBd09(p.getLon(), p.getLat()))
                .collect(Collectors.toList());
    }

    /**
     * 批量将BD09坐标点列表直接转换为WGS84坐标点列表
     */
    public List<CoordinatePoint> bd09ToWgs84(List<CoordinatePoint> points) {
        return points.stream()
                .map(p -> bd09ToWgs84(p.getLon(), p.getLat()))
                .collect(Collectors.toList());
    }

    /**
     * 批量将TrackPoint列表从WGS84坐标转换为GCJ02坐标（修改原对象）
     */
    public List<TrackPoint> wgs84ToGcj02TrackPoints(List<TrackPoint> points) {
        return points.stream()
                .peek(p -> {
                    CoordinatePoint result = wgs84ToGcj02(p.getLon(), p.getLat());
                    p.setLon(result.getLon());
                    p.setLat(result.getLat());
                })
                .collect(Collectors.toList());
    }

    /**
     * 批量将TrackPoint列表从GCJ02坐标转换为WGS84坐标（修改原对象）
     */
    public List<TrackPoint> gcj02ToWgs84TrackPoints(List<TrackPoint> points) {
        return points.stream()
                .peek(p -> {
                    CoordinatePoint result = gcj02ToWgs84(p.getLon(), p.getLat());
                    p.setLon(result.getLon());
                    p.setLat(result.getLat());
                })
                .collect(Collectors.toList());
    }

    /**
     * 批量将TrackPoint列表从GCJ02坐标转换为BD09坐标（创建新对象）
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
     */
    public boolean outOfChina(double lon, double lat) {
        return lon < 72.004 || lon > 137.8347 || lat < 0.8293 || lat > 55.8271;
    }

}
