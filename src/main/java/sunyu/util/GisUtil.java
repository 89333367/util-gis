package sunyu.util;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;
import java.util.stream.Collectors;

import org.geotools.geometry.jts.JTS;
import org.geotools.referencing.CRS;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Envelope;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.GeometryCollection;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.LineString;
import org.locationtech.jts.geom.MultiPolygon;
import org.locationtech.jts.geom.Polygon;
import org.locationtech.jts.geom.prep.PreparedGeometry;
import org.locationtech.jts.geom.prep.PreparedGeometryFactory;
import org.locationtech.jts.io.WKTReader;
import org.locationtech.jts.operation.buffer.BufferOp;
import org.locationtech.jts.operation.buffer.BufferParameters;
import org.locationtech.jts.operation.union.UnaryUnionOp;
import org.locationtech.jts.simplify.DouglasPeuckerSimplifier;
import org.opengis.referencing.crs.CoordinateReferenceSystem;
import org.opengis.referencing.operation.MathTransform;

import cn.hutool.core.bean.BeanUtil;
import cn.hutool.log.Log;
import cn.hutool.log.LogFactory;
import sunyu.util.pojo.CoordinatePoint;
import sunyu.util.pojo.OutlinePart;
import sunyu.util.pojo.SplitRoadResult;
import sunyu.util.pojo.TrackPoint;

/**
 * GIS工具类
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
        log.info("[构建{}] 开始", this.getClass().getSimpleName());
        // 其他初始化语句（预留扩展点）
        // 记录工具类构建结束日志
        log.info("[构建{}] 结束", this.getClass().getSimpleName());
        // 保存配置参数引用
        this.config = config;
    }

    /**
     * 内部配置类，定义一些常量
     * 使用内部类封装配置参数，提高代码组织性和封装性
     */
    private static class Config {
        // 最小亩数阈值，过滤小面积多边形
        private final double MIN_MU_THRESHOLD = 0.7;

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

        // 默认轮廓返回的最多多边形数量（TopN）
        private final int DEFAULT_MAX_OUTLINE_SEGMENTS = 10;

        private final int DEFAULT_UNION_GROUP_SIZE = 600; // 分组合并组大小：600，避免 GeometryCollection 一次性过大
        private final int DEFAULT_BUFFER_QUADRANT = 4; // 圆近似细分：4，提升性能，误差满足道路宽度场景
        private final int DEFAULT_BUCKET_TARGET = 300; // 目标桶数：300，平衡桶内复杂度与最终合并规模
        private final int DEFAULT_BUCKET_CELL_MIN_FACTOR = 5; // cellSize 下限系数：5，避免过多小桶导致合并碎片化
        private final int DEFAULT_BUCKET_CELL_MAX_FACTOR = 36; // cellSize 上限系数：36

        // 线缓冲简化容差参数（基于宽度因子与最小值）
        private final double LINE_SIMPLIFY_FACTOR = 0.25; // 宽度因子，越小越保留细节
        private final double LINE_SIMPLIFY_MIN = 0.2; // 最小容差（米）

        private volatile MathTransform txWgsToMercator;
        private volatile MathTransform txMercatorToWgs;
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
     * 几何图形合并工具（内部类）
     * 使用分组批量合并策略避免大量单独 union 操作的性能问题
     */
    private Geometry unionGeometries(List<Geometry> geometries) {
        if (geometries == null || geometries.isEmpty()) {
            return null;
        }
        if (geometries.size() == 1) {
            return geometries.get(0);
        }

        // 优先使用 CascadedPolygonUnion 合并多边形，性能更优
        java.util.List<Polygon> polys = new java.util.ArrayList<>();
        for (Geometry g : geometries) {
            if (g instanceof Polygon) {
                polys.add((Polygon) g);
            } else if (g instanceof MultiPolygon) {
                MultiPolygon mp = (MultiPolygon) g;
                for (int i = 0; i < mp.getNumGeometries(); i++) {
                    polys.add((Polygon) mp.getGeometryN(i));
                }
            }
        }
        if (!polys.isEmpty()) {
            return UnaryUnionOp.union(polys);
        }

        // 非多边形集合，退回 GeometryCollection.union() 或分组合并
        int groupSize = config.DEFAULT_UNION_GROUP_SIZE;
        if (geometries.size() <= groupSize) {
            GeometryCollection collection = config.geometryFactory.createGeometryCollection(
                    geometries.toArray(new Geometry[0]));
            return collection.union();
        }
        return unionInGroups(geometries, groupSize);
    }

    private Geometry unionInGroups(List<Geometry> geometries, int groupSize) {
        if (geometries.size() <= groupSize) {
            GeometryCollection collection = config.geometryFactory.createGeometryCollection(
                    geometries.toArray(new Geometry[0]));
            return collection.union();
        }
        List<Geometry> groups = new java.util.ArrayList<>();
        for (int i = 0; i < geometries.size(); i += groupSize) {
            int end = Math.min(i + groupSize, geometries.size());
            List<Geometry> group = geometries.subList(i, end);
            GeometryCollection collection = config.geometryFactory.createGeometryCollection(
                    group.toArray(new Geometry[0]));
            groups.add(collection.union());
        }
        return unionInGroups(groups, groupSize);
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
                throw new RuntimeException("解析CRS出错: " + key, e);
            }
        });
    }

    private MathTransform getTxWgsToMercator() {
        MathTransform tx = config.txWgsToMercator;
        if (tx == null) {
            synchronized (this) {
                if (config.txWgsToMercator == null) {
                    try {
                        CoordinateReferenceSystem src = getCachedCRS(config.WGS84);
                        CoordinateReferenceSystem tgt = getCachedCRS(config.WEB_MERCATOR);
                        config.txWgsToMercator = CRS.findMathTransform(src, tgt, true);
                    } catch (Exception e) {
                        throw new RuntimeException("初始化WGS->Mercator变换失败", e);
                    }
                }
                tx = config.txWgsToMercator;
            }
        }
        return tx;
    }

    private MathTransform getTxMercatorToWgs() {
        MathTransform tx = config.txMercatorToWgs;
        if (tx == null) {
            synchronized (this) {
                if (config.txMercatorToWgs == null) {
                    try {
                        CoordinateReferenceSystem src = getCachedCRS(config.WEB_MERCATOR);
                        CoordinateReferenceSystem tgt = getCachedCRS(config.WGS84);
                        config.txMercatorToWgs = CRS.findMathTransform(src, tgt, true);
                    } catch (Exception e) {
                        throw new RuntimeException("初始化Mercator->WGS变换失败", e);
                    }
                }
                tx = config.txMercatorToWgs;
            }
        }
        return tx;
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
        // 通过 getTxWgsToMercator() 懒加载，无需重复初始化逻辑
        return JTS.transform(g, getTxWgsToMercator());
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
            // 通过 getTxMercatorToWgs() 懒加载，无需重复初始化逻辑
            Geometry result = JTS.transform(g, getTxMercatorToWgs());
            if (!isValidWgs84Geometry(result)) {
                CoordinateReferenceSystem src = getCachedCRS(config.WEB_MERCATOR);
                CoordinateReferenceSystem tgt = getCachedCRS(config.WGS84);
                MathTransform strict = CRS.findMathTransform(src, tgt, false);
                result = JTS.transform(g, strict);
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

    // 统一过滤与排序 TrackPoint
    private List<TrackPoint> filterAndSortTrackPoints(List<TrackPoint> seg, boolean requireTime) {
        if (seg == null)
            return java.util.Collections.emptyList();
        java.util.Comparator<java.time.LocalDateTime> cmp = java.util.Comparator
                .nullsLast(java.util.Comparator.naturalOrder());
        return seg.stream()
                .filter(p -> !requireTime || p.getTime() != null)
                .filter(p -> Math.abs(p.getLon()) <= 180 && Math.abs(p.getLat()) <= 90)
                .filter(p -> !(p.getLon() == 0 && p.getLat() == 0))
                .sorted(java.util.Comparator.comparing(TrackPoint::getTime, cmp))
                .collect(java.util.stream.Collectors.toList());
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
    private Geometry buildOutlineCore(List<TrackPoint> seg, double totalWidthM) throws Exception {
        long startTime = System.currentTimeMillis();
        double widthM = totalWidthM / 2.0;
        log.trace("[buildOutlineCore] 开始处理轨迹轮廓，轨迹点数: {}, 宽度(单侧): {} 米",
                seg != null ? seg.size() : 0, widthM);

        if (seg == null) {
            throw new IllegalArgumentException("轨迹段至少需要3个点");
        }
        if (widthM < 0) {
            throw new IllegalArgumentException("所有参数必须为非负数");
        }

        long processStartTime = System.currentTimeMillis();
        List<TrackPoint> sortedSeg = filterAndSortTrackPoints(seg, false);
        long processTime = System.currentTimeMillis() - processStartTime;
        log.trace("[buildOutlineCore] 轨迹点过滤完成，原始点数: {}, 过滤后点数: {}, 过滤耗时: {}ms",
                seg.size(), sortedSeg.size(), processTime);

        if (sortedSeg.size() < 3) {
            throw new IllegalArgumentException("轨迹段至少需要3个有效点");
        }

        try {
            double minLon = sortedSeg.stream().mapToDouble(TrackPoint::getLon).min().orElse(0);
            double maxLon = sortedSeg.stream().mapToDouble(TrackPoint::getLon).max().orElse(0);
            double minLat = sortedSeg.stream().mapToDouble(TrackPoint::getLat).min().orElse(0);
            double maxLat = sortedSeg.stream().mapToDouble(TrackPoint::getLat).max().orElse(0);

            log.trace("[buildOutlineCore] 轨迹范围: 经度[{}, {}], 纬度[{}, {}]", minLon, maxLon, minLat, maxLat);

            long buildStartTime = System.currentTimeMillis();
            Geometry result;
            if (sortedSeg.size() <= 200) {
                log.trace("[buildOutlineCore] 选择LineBuffer策略（点数:{}）", sortedSeg.size());
                result = buildOutlineByLineBuffer(sortedSeg, widthM);
            } else {
                log.trace("[buildOutlineCore] 选择PointBuffer策略（点数:{}）", sortedSeg.size());
                result = buildOutlineBySimpleBuffers(sortedSeg, widthM);
            }
            long buildTime = System.currentTimeMillis() - buildStartTime;

            // 基于宽度的最小面积阈值（直接使用宽度，不再依赖配置因子）
            double minAreaThresholdM2 = Math.PI * widthM * widthM;
            if (result instanceof MultiPolygon) {
                MultiPolygon mp = (MultiPolygon) result;
                List<Polygon> keep = new ArrayList<>();
                for (int i = 0; i < mp.getNumGeometries(); i++) {
                    Polygon poly = (Polygon) mp.getGeometryN(i);
                    Geometry projPoly = wgs84ToWebMercator(poly);
                    double areaM2 = projPoly.getArea();
                    if (areaM2 >= minAreaThresholdM2) {
                        keep.add(poly);
                    }
                }
                if (!keep.isEmpty()) {
                    result = config.geometryFactory.createMultiPolygon(keep.toArray(new Polygon[0]));
                } else {
                    log.trace("[buildOutlineCore] 过滤后无有效多边形，阈值面积: {} m²", minAreaThresholdM2);
                    result = config.geometryFactory.createMultiPolygon(new Polygon[0]);
                }
                log.trace("[buildOutlineCore] 小多边形过滤完成：原 {} 个 → 保留 {} 个；阈值面积: {} m²",
                        mp.getNumGeometries(), keep.size(), minAreaThresholdM2);
            }

            log.trace("[buildOutlineCore] 轮廓构建完成，构建耗时: {}ms", buildTime);

            if (result instanceof MultiPolygon && result.getNumGeometries() > 1) {
                log.debug("[buildOutlineCore] 输入参数 点数={} 总宽度={}m 单侧宽度={}m", seg != null ? seg.size() : 0, totalWidthM,
                        widthM);
                log.debug("[buildOutlineCore] 结果为MultiPolygon，包含 {} 个部分", result.getNumGeometries());
                long endTime = System.currentTimeMillis();
                log.trace("[buildOutlineCore] 总计耗时: {}ms", endTime - startTime);
                return result;
            }

            if (result instanceof Polygon) {
                log.trace("[buildOutlineCore] 结果为Polygon");
                long endTime = System.currentTimeMillis();
                log.trace("[buildOutlineCore] 总计耗时: {}ms", endTime - startTime);
                return result;
            }

            if (result instanceof MultiPolygon && result.getNumGeometries() == 1) {
                result = (Polygon) ((MultiPolygon) result).getGeometryN(0);
                log.trace("[buildOutlineCore] 结果为单个Polygon（由单个MultiPolygon扁平化）");
                long endTime = System.currentTimeMillis();
                log.trace("[buildOutlineCore] 总计耗时: {}ms", endTime - startTime);
                return result;
            }

            long endTime = System.currentTimeMillis();
            log.info("[buildOutlineCore] 处理完成，总计耗时: {}ms", endTime - startTime);

            return result;
        } catch (Exception e) {
            log.error("构建轨迹轮廓失败: " + e.getMessage(), e);
            throw new Exception("构建轨迹轮廓失败: " + e.getMessage(), e);
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
        log.trace("[buildOutlineBySimpleBuffers] 开始处理 {} 个轨迹点，缓冲区宽度: {} 米", seg.size(), widthM);

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

            // 创建缓冲区（降低圆近似复杂度以提升性能）
            long startBuffer = System.currentTimeMillis();
            BufferParameters params = new BufferParameters();
            int quadSeg = config.DEFAULT_BUFFER_QUADRANT;
            params.setQuadrantSegments(quadSeg); // 使用配置中的圆弧细分段数
            // 点缓冲固定使用圆端+圆角，避免意外退化
            params.setEndCapStyle(BufferParameters.CAP_ROUND);
            params.setJoinStyle(BufferParameters.JOIN_ROUND);
            Geometry buffer = BufferOp.bufferOp(projPoint, widthM, params);
            bufferTime += System.currentTimeMillis() - startBuffer;

            pointBuffers.add(buffer);

            totalTime += System.currentTimeMillis() - startPoint;

            // 每处理100个点打印一次进度
            if (i > 0 && i % 100 == 0) {
                log.trace("[buildOutlineBySimpleBuffers] 已处理 {}/{} 个点, 转换耗时: {}ms, 缓冲耗时: {}ms, 平均每个点耗时: {}ms",
                        i, seg.size(), convertTime, bufferTime, totalTime / i);
            }
        }

        log.trace("[buildOutlineBySimpleBuffers] 点缓冲区创建完成，共 {} 个缓冲区, 转换总耗时: {}ms, 缓冲总耗时: {}ms",
                pointBuffers.size(), convertTime, bufferTime);

        // 使用GeometryCollection优化合并所有缓冲区
        log.trace("[buildOutlineBySimpleBuffers] 开始使用GeometryCollection优化合并 {} 个缓冲区", pointBuffers.size());

        // 分桶：在 WebMercator 上按网格将缓冲分组，降低相互参与的多边形数量
        long startBucket = System.currentTimeMillis();
        // 基于整体包络自适应格子大小，目标桶数约300
        Envelope overall = new Envelope();
        for (Geometry g : pointBuffers) {
            overall.expandToInclude(g.getEnvelopeInternal());
        }
        double area = overall.getWidth() * overall.getHeight();
        double cellSize;
        if (area > 0) {
            int targetBuckets = config.DEFAULT_BUCKET_TARGET;
            int minFactor = config.DEFAULT_BUCKET_CELL_MIN_FACTOR;
            int maxFactor = config.DEFAULT_BUCKET_CELL_MAX_FACTOR;
            targetBuckets = Math.max(1, targetBuckets);
            minFactor = Math.max(1, minFactor);
            maxFactor = Math.max(minFactor, maxFactor);
            cellSize = Math.sqrt(area / targetBuckets);
            cellSize = Math.max(widthM * minFactor, Math.min(cellSize, widthM * maxFactor));
            log.trace("[buildOutlineBySimpleBuffers] 分桶参数 targetBuckets={} minFactor={} maxFactor={} cellSize={}",
                    targetBuckets, minFactor, maxFactor, Math.round(cellSize));
        } else {
            cellSize = widthM * 16;
        }
        Map<String, List<Geometry>> buckets = new HashMap<>();
        for (Geometry g : pointBuffers) {
            Envelope e = g.getEnvelopeInternal();
            double cx = e.getMinX() + e.getWidth() / 2.0;
            double cy = e.getMinY() + e.getHeight() / 2.0;
            long bx = (long) Math.floor(cx / cellSize);
            long by = (long) Math.floor(cy / cellSize);
            String key = bx + ":" + by;
            buckets.computeIfAbsent(key, k -> new ArrayList<>()).add(g);
        }
        long bucketTime = System.currentTimeMillis() - startBucket;
        int bucketCount = buckets.size();
        int avgBucketSize = bucketCount == 0 ? 0 : (pointBuffers.size() / bucketCount);
        log.trace("[buildOutlineBySimpleBuffers] 分桶完成 桶数={} 平均每桶={} 耗时={}ms", bucketCount, avgBucketSize, bucketTime);

        // 桶内合并（并行）
        long startBucketUnion = System.currentTimeMillis();
        List<Geometry> bucketUnions = buckets.entrySet()
                .parallelStream()
                .map(e -> unionGeometries(e.getValue()))
                .filter(java.util.Objects::nonNull)
                .collect(java.util.stream.Collectors.toList());
        long bucketUnionTime = System.currentTimeMillis() - startBucketUnion;
        log.trace("[buildOutlineBySimpleBuffers] 桶内合并完成(并行) 桶数={} 合并结果数={} 耗时={}ms", bucketCount, bucketUnions.size(),
                bucketUnionTime);

        // 最终合并（合并各桶的结果）
        long startFinalUnion = System.currentTimeMillis();
        Geometry union;
        if (bucketUnions.isEmpty()) {
            GeometryCollection collection = config.geometryFactory
                    .createGeometryCollection(pointBuffers.toArray(new Geometry[0]));
            union = collection.union();
        } else {
            union = unionGeometries(bucketUnions);
        }
        long finalUnionTime = System.currentTimeMillis() - startFinalUnion;
        long unionTime = bucketUnionTime + finalUnionTime;
        log.trace("[buildOutlineBySimpleBuffers] 最终合并完成 桶内:{}ms 总合并:{}ms", bucketUnionTime, finalUnionTime);

        // 转换回WGS84坐标系
        long startBackConvert = System.currentTimeMillis();
        Geometry result = webMercatorToWgs84(union);
        long backConvertTime = System.currentTimeMillis() - startBackConvert;

        long endTime = System.currentTimeMillis();
        log.trace("[buildOutlineBySimpleBuffers] 处理完成，总计耗时: {}ms (转换:{}ms, 缓冲:{}ms, 合并:{}ms, 回转:{}ms)",
                endTime - startTime, convertTime, bufferTime, unionTime, backConvertTime);

        return result;
    }

    private Geometry buildOutlineByLineBuffer(List<TrackPoint> seg, double widthM) throws Exception {
        long startTime = System.currentTimeMillis();
        log.trace("[buildOutlineByLineBuffer] 开始处理 {} 个轨迹点，缓冲区宽度: {} 米", seg.size(), widthM);

        Coordinate[] coords = new Coordinate[seg.size()];
        for (int i = 0; i < seg.size(); i++) {
            TrackPoint p = seg.get(i);
            coords[i] = new Coordinate(p.getLon(), p.getLat());
        }
        LineString line = config.geometryFactory.createLineString(coords);

        long t1 = System.currentTimeMillis();
        Geometry projLine = wgs84ToWebMercator(line);
        long convertTime = System.currentTimeMillis() - t1;

        // 简化容差可配置，避免过度圆滑导致桥接
        double factor = config.LINE_SIMPLIFY_FACTOR;
        double minTol = config.LINE_SIMPLIFY_MIN;
        double tolerance = Math.max(minTol, widthM * factor);
        long t2 = System.currentTimeMillis();
        Geometry simpleLine = DouglasPeuckerSimplifier.simplify(projLine, tolerance);
        long simplifyTime = System.currentTimeMillis() - t2;

        long t3 = System.currentTimeMillis();
        // 线缓冲使用可配置的端帽/连接样式，减少拐点外扩
        BufferParameters lp = new BufferParameters();
        int quadSeg2 = config.DEFAULT_BUFFER_QUADRANT;
        lp.setQuadrantSegments(quadSeg2);
        lp.setEndCapStyle(BufferParameters.CAP_FLAT);
        lp.setJoinStyle(BufferParameters.JOIN_BEVEL);
        lp.setMitreLimit(2.0);
        Geometry buffered = BufferOp.bufferOp(simpleLine, widthM, lp);
        long bufferTime = System.currentTimeMillis() - t3;

        long t4 = System.currentTimeMillis();
        Geometry result = webMercatorToWgs84(buffered);
        long backConvertTime = System.currentTimeMillis() - t4;

        long endTime = System.currentTimeMillis();
        log.trace("[buildOutlineByLineBuffer] 处理完成，总计耗时: {}ms (投影:{}ms, 简化:{}ms, 缓冲:{}ms, 回转:{}ms)",
                endTime - startTime, convertTime, simplifyTime, bufferTime, backConvertTime);

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
        dLat = (dLat * 180.0)
                / ((config.semiMajorAxis * (1 - config.eccentricitySquared)) / (magic * sqrtMagic) * Math.PI);
        dLon = (dLon * 180.0) / (config.semiMajorAxis / sqrtMagic * Math.cos(radLat) * Math.PI);
        return new double[] { dLon, dLat };
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

    // 新增：按面积倒序保留最大的 N 个多边形
    private Geometry keepLargestPolygons(Geometry geometry, int maxCount) {
        if (geometry == null)
            return null;
        if (maxCount <= 0)
            return geometry;

        GeometryFactory gf = geometry.getFactory();
        java.util.List<Polygon> polys = new java.util.ArrayList<>();

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
        polys.sort(java.util.Comparator.comparingDouble(Polygon::getArea).reversed());

        int limit = Math.min(maxCount, polys.size());
        if (limit == 1) {
            return polys.get(0);
        }
        Polygon[] top = polys.subList(0, limit).toArray(new Polygon[0]);
        return gf.createMultiPolygon(top);
    }

    // 新增：按最小亩数阈值过滤小面积多边形
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
            java.util.List<Polygon> kept = new java.util.ArrayList<>();
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

    // 球面环面积（平方米），与 Turf.js 的 ringArea 保持一致
    private static double ringArea(org.locationtech.jts.geom.LineString ring) {
        org.locationtech.jts.geom.Coordinate[] coords = ring.getCoordinates();
        int len = (coords == null) ? 0 : coords.length;
        if (len <= 2)
            return 0.0;
        double area = 0.0;
        for (int i = 0; i < len; i++) {
            org.locationtech.jts.geom.Coordinate p1, p2, p3;
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

    private static double toRad(double deg) {
        return deg * Math.PI / 180.0;
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
     * 计算几何图形的面积（mu单位，支持 Polygon 与 MultiPolygon），按 WGS84 球面公式对齐 Turf.js
     *
     * @param outline 几何图形（POLYGON 或 MULTIPOLYGON）
     *
     * @return 面积（mu），以亩为单位，保留4位小数
     *
     * @throws RuntimeException 如果面积计算过程中发生错误
     */
    public double calcMu(Geometry outline) throws RuntimeException {
        try {
            double areaSqm = 0.0;
            if (outline instanceof org.locationtech.jts.geom.Polygon) {
                org.locationtech.jts.geom.Polygon p = (org.locationtech.jts.geom.Polygon) outline;
                double areaOuter = Math.abs(ringArea(p.getExteriorRing()));
                double holesArea = 0.0;
                for (int i = 0; i < p.getNumInteriorRing(); i++) {
                    holesArea += Math.abs(ringArea(p.getInteriorRingN(i)));
                }
                areaSqm = areaOuter - holesArea;
            } else if (outline instanceof org.locationtech.jts.geom.MultiPolygon) {
                org.locationtech.jts.geom.MultiPolygon mp = (org.locationtech.jts.geom.MultiPolygon) outline;
                for (int i = 0; i < mp.getNumGeometries(); i++) {
                    org.locationtech.jts.geom.Polygon p = (org.locationtech.jts.geom.Polygon) mp.getGeometryN(i);
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
     * WKT计算亩数（WGS84坐标系），支持POLYGON与MULTIPOLYGON
     * 直接计算WGS84坐标系下几何图形的面积，结果以亩为单位
     *
     * @param wkt WKT字符串，要求为WGS84坐标系
     *
     * @return 面积（mu），以亩为单位，保留4位小数
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
            Geometry polygon = fromWkt(polygonWkt);
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
     * 判断两个多边形是否相交（支持POLYGON与MULTIPOLYGON）
     *
     * @param wktA 第一个几何的WKT（WGS84）
     * @param wktB 第二个几何的WKT（WGS84）
     * @return 是否相交
     * @throws Exception 解析或计算过程中出错
     */
    public boolean intersects(String wktA, String wktB) throws Exception {
        Geometry g1 = fromWkt(wktA);
        Geometry g2 = fromWkt(wktB);
        return g1.intersects(g2);
    }

    /**
     * 判断两个多边形是否相邻（共享边界，支持POLYGON与MULTIPOLYGON）
     *
     * @param wktA 第一个几何的WKT（WGS84）
     * @param wktB 第二个几何的WKT（WGS84）
     * @return 是否相邻（touches）
     * @throws Exception 解析或计算过程中出错
     */
    public boolean touches(String wktA, String wktB) throws Exception {
        Geometry g1 = fromWkt(wktA);
        Geometry g2 = fromWkt(wktB);
        return g1.touches(g2);
    }

    /**
     * 创建缓冲区
     */
    public Geometry buffer(Geometry geom, double distance, double originLon) throws Exception {
        try {
            // 转换到Web Mercator投影坐标
            Geometry proj = wgs84ToWebMercator(geom);
            Geometry buffer = proj.buffer(distance);
            // 转换回WGS84坐标系
            return webMercatorToWgs84(buffer);
        } catch (Exception e) {
            throw new Exception("创建缓冲区时出错: " + e.getMessage(), e);
        }
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
                    TrackPoint tp = new TrackPoint(p.getTime(), result.getLon(), result.getLat());
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
                    CoordinatePoint result = bd09ToGcj02(p.getLon(), p.getLat());
                    TrackPoint tp = new TrackPoint(p.getTime(), result.getLon(), result.getLat());
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
                    CoordinatePoint result = wgs84ToBd09(p.getLon(), p.getLat());
                    TrackPoint tp = new TrackPoint(p.getTime(), result.getLon(), result.getLat());
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
                    CoordinatePoint result = bd09ToWgs84(p.getLon(), p.getLat());
                    TrackPoint tp = new TrackPoint(p.getTime(), result.getLon(), result.getLat());
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

    /**
     * 基于轨迹点与总宽度生成轮廓，并按面积保留 Top-N 最大区块。
     *
     * <p>参数说明</p>
     * - seg：轨迹点列表（WGS84，经度范围[-180,180]、纬度范围[-90,90]），至少 3 个有效点；
     *   方法内部会按时间升序处理，并自动过滤非法点（越界以及 0,0）。
     * - totalWidthM：道路“总宽度”（米），等于左右两侧宽度之和；内部以 totalWidthM/2 作为单侧缓冲宽度。
     *
     * <p>返回结果</p>
     * - 仅形成一个区块时返回 `Polygon`；形成多个区块时返回 `MultiPolygon`；
     * - 会按面积从大到小保留前 N 个区块（默认 N=10）；
     * - 若过滤后无有效区块，则返回空 `MultiPolygon`；
     * - 返回坐标系为 WGS84。
     *
     * <p>行为特性</p>
     * - 以宽度为依据的最小面积阈值过滤小区块，阈值约为 π*(totalWidthM/2)^2；
     * - 面积倒序裁剪，最多保留 `maxSegments` 个区块；若仅 1 个区块将返回 `Polygon`。
     *
     * <p>异常</p>
     * - IllegalArgumentException：seg 少于 3 个有效点，或 totalWidthM < 0；
     * - Exception：坐标转换或几何运算异常。
     *
     * <p>示例</p>
     * `GisUtil util = GisUtil.builder().build();`
     * `Geometry outline = util.splitRoad(points, 8.0);`
     *
     * @param seg         轨迹点列表（WGS84），至少 3 个有效点
     * @param totalWidthM 道路总宽度（米），左右合计，必须非负
     * @return            轮廓几何，可能为 `Polygon` 或 `MultiPolygon`（裁剪后最多保留前 N 个）
     * @throws Exception  坐标转换或几何运算异常
     */
    public SplitRoadResult splitRoad(List<TrackPoint> seg, double totalWidthM) throws Exception {
        return splitRoad(seg, totalWidthM, config.DEFAULT_MAX_OUTLINE_SEGMENTS);
    }

    /**
     * 基于轨迹点与总宽度生成轮廓，并按面积保留不超过 maxSegments 的最大区块。
     *
     * <p>参数说明</p>
     * - seg：轨迹点列表（WGS84，经度范围[-180,180]、纬度范围[-90,90]），至少 3 个有效点；
     *   方法内部会按时间升序处理，并自动过滤非法点（越界以及 0,0）。
     * - totalWidthM：道路“总宽度”（米），等于左右两侧宽度之和；内部以 totalWidthM/2 作为单侧缓冲宽度。
     * - maxSegments：返回区块数量上限；为 null 或 <=0 时使用默认值 N=10；
     *   当原始区块数 ≤ 上限时不会进行裁剪（不输出“裁剪结果”日志）。
     *
     * <p>返回结果</p>
     * - 仅形成一个区块时返回 `Polygon`；形成多个区块时返回 `MultiPolygon`；
     * - 若过滤后无有效区块，则返回空 `MultiPolygon`；
     * - 返回坐标系为 WGS84。
     *
     * <p>行为特性</p>
     * - 以宽度为依据的最小面积阈值过滤小区块，阈值约为 π*(totalWidthM/2)^2；
     * - 面积倒序裁剪，最多保留 `maxSegments` 个区块；若仅 1 个区块将返回 `Polygon`。
     *
     * <p>异常</p>
     * - IllegalArgumentException：seg 少于 3 个有效点，或 totalWidthM < 0；
     * - Exception：坐标转换或几何运算异常。
     *
     * <p>示例</p>
     * `GisUtil util = GisUtil.builder().build();`
     * `Geometry outline = util.splitRoad(points, 8.0, 5);`
     *
     * @param seg         轨迹点列表（WGS84），至少 3 个有效点
     * @param totalWidthM 道路总宽度（米），左右合计，必须非负
     * @param maxSegments 返回的区块数量上限；为 null 或 <=0 时取默认值 N=10
     * @return            轮廓几何，可能为 `Polygon` 或 `MultiPolygon`（最多保留 `maxSegments` 个）
     * @throws Exception  坐标转换或几何运算异常
     */
    public SplitRoadResult splitRoad(List<TrackPoint> seg, double totalWidthM, Integer maxSegments) throws Exception {
        long t0 = System.currentTimeMillis();
        Geometry outline = buildOutlineCore(seg, totalWidthM);
        long t1 = System.currentTimeMillis();
        log.debug("[splitRoad] outline构建完成 type={} parts={}", outline.getGeometryType(),
                (outline instanceof MultiPolygon) ? outline.getNumGeometries() : 1);
        log.debug("[splitRoad] outline构建耗时: {}ms", (t1 - t0));
        int partsOutline = (outline instanceof MultiPolygon) ? outline.getNumGeometries() : 1;

        int limit = (maxSegments == null || maxSegments <= 0) ? config.DEFAULT_MAX_OUTLINE_SEGMENTS
                : maxSegments.intValue();
        log.debug("[splitRoad] 输入参数 点数={} 总宽度={}m 返回上限={}", seg != null ? seg.size() : 0, totalWidthM, limit);
        long tTrimStart = System.currentTimeMillis();
        Geometry trimmed = keepLargestPolygons(outline, limit);
        int partsAfterKeep = (trimmed instanceof MultiPolygon) ? trimmed.getNumGeometries() : 1;
        // 依据最小亩数阈值过滤小区块
        trimmed = removeSmallMuPolygons(trimmed, config.MIN_MU_THRESHOLD);
        int partsAfterFilter = (trimmed instanceof MultiPolygon) ? trimmed.getNumGeometries() : 1;
        int removedByArea = partsAfterKeep - partsAfterFilter;
        long tTrimEnd = System.currentTimeMillis();
        log.trace("[splitRoad] 裁剪+小面积过滤完成 type={} parts={} 移除={} 阈值={}mu 耗时={}ms", trimmed.getGeometryType(),
                partsAfterFilter, removedByArea, config.MIN_MU_THRESHOLD, (tTrimEnd - tTrimStart));

        // 最终返回日志：当少于上限时说明原因，否则保持原样
        if (partsAfterFilter < limit) {
            String reason;
            if (partsAfterFilter != partsAfterKeep) {
                reason = "小面积过滤";
            } else if (partsOutline < limit) {
                reason = "原始区块数小于上限";
            } else {
                reason = "面积裁剪后区块数不超过上限";
            }
            log.debug("[splitRoad] 返回结果 type={} parts={} (上限={}，原因：{}；小面积移除={})", trimmed.getGeometryType(),
                    partsAfterFilter, limit, reason, removedByArea);
        } else {
            log.debug("[splitRoad] 返回结果 type={} parts={}", trimmed.getGeometryType(), partsAfterFilter);
        }

        // 准备有效的轨迹点（与轮廓构建一致的过滤逻辑）
        long tFilterStart = System.currentTimeMillis();
        List<TrackPoint> sortedSeg = filterAndSortTrackPoints(seg, true);
        long tFilterEnd = System.currentTimeMillis();
        log.trace("[splitRoad] 点过滤+排序完成 有效点={} 耗时={}ms", sortedSeg.size(), (tFilterEnd - tFilterStart));

        java.util.List<OutlinePart> parts = new java.util.ArrayList<>();

        if (trimmed instanceof Polygon) {
            log.trace("[splitRoad] 进入Polygon分支");
            Polygon poly = (Polygon) trimmed;
            PreparedGeometry preparedPoly = PreparedGeometryFactory.prepare(poly);
            Envelope env = poly.getEnvelopeInternal();
            // 找出属于该区块的点（点在多边形内）
            long tWithinStart = System.currentTimeMillis();
            java.util.List<TrackPoint> inPoly = new java.util.ArrayList<>();
            for (TrackPoint p : sortedSeg) {
                Coordinate c = new Coordinate(p.getLon(), p.getLat());
                if (!env.contains(c)) {
                    continue;
                }
                Geometry point = config.geometryFactory.createPoint(c);
                if (preparedPoly.contains(point)) {
                    inPoly.add(p);
                }
            }
            long tWithinEnd = System.currentTimeMillis();
            log.trace("[splitRoad] Polygon点过滤完成 inPoly={} 耗时={}ms", inPoly.size(), (tWithinEnd - tWithinStart));

            java.time.LocalDateTime start = inPoly.isEmpty() ? null : inPoly.get(0).getTime();
            java.time.LocalDateTime end = inPoly.isEmpty() ? null : inPoly.get(inPoly.size() - 1).getTime();
            log.trace("[splitRoad] Polygon 时间范围 start={} end={} inPoly={}", start, end, inPoly.size());

            long tMuStart = System.currentTimeMillis();
            double mu = calcMu(poly);
            long tMuEnd = System.currentTimeMillis();
            log.trace("[splitRoad] Polygon mu计算完成 值={} 耗时={}ms", mu, (tMuEnd - tMuStart));

            long tWktStart = System.currentTimeMillis();
            String wkt = toWkt(poly);
            long tWktEnd = System.currentTimeMillis();
            log.trace("[splitRoad] Polygon WKT生成完成 长度={} 耗时={}ms", wkt.length(), (tWktEnd - tWktStart));

            parts.add(new OutlinePart(poly, start, end, mu, wkt));
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
                java.util.List<java.util.List<TrackPoint>> sessions = new java.util.ArrayList<>();
                java.util.List<TrackPoint> current = null;
                for (TrackPoint p : sortedSeg) {
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
                            current = new java.util.ArrayList<>();
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

                java.time.LocalDateTime start = null;
                java.time.LocalDateTime end = null;
                if (!sessions.isEmpty()) {
                    // 使用进入该区块的最早时间和离开该区块的最晚时间
                    sessions.sort(java.util.Comparator.comparing(
                            s -> s.get(0).getTime(),
                            java.util.Comparator.nullsLast(java.util.Comparator.naturalOrder())));
                    java.util.List<TrackPoint> firstSess = sessions.get(0);
                    java.util.List<TrackPoint> lastSess = sessions.get(sessions.size() - 1);
                    start = firstSess.isEmpty() ? null : firstSess.get(0).getTime();
                    end = lastSess.isEmpty() ? null : lastSess.get(lastSess.size() - 1).getTime();
                }
                log.trace("[splitRoad] Part#{} 时间范围 start={} end={} (sessions={})", i, start, end, sessions.size());

                long tMuStart = System.currentTimeMillis();
                double mu = calcMu(poly);
                long tMuEnd = System.currentTimeMillis();
                log.trace("[splitRoad] Part#{} mu计算完成 值={} 耗时={}ms", i, mu, (tMuEnd - tMuStart));

                long tWktStart = System.currentTimeMillis();
                String wkt = toWkt(poly);
                long tWktEnd = System.currentTimeMillis();
                log.trace("[splitRoad] Part#{} WKT生成完成 长度={} 耗时={}ms", i, wkt.length(), (tWktEnd - tWktStart));

                parts.add(new OutlinePart(poly, start, end, mu, wkt));
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
        log.debug("[splitRoad] 总耗时={}ms", (tTotalEnd - t0));
        return new SplitRoadResult(trimmed, parts, outlineWkt);
    }

}
