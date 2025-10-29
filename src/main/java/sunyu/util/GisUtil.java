package sunyu.util;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;

import org.geotools.geometry.jts.JTS;
import org.geotools.referencing.CRS;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Envelope;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.GeometryCollection;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.MultiPolygon;
import org.locationtech.jts.geom.Polygon;
import org.locationtech.jts.geom.prep.PreparedGeometry;
import org.locationtech.jts.geom.prep.PreparedGeometryFactory;
import org.locationtech.jts.io.WKTReader;
import org.locationtech.jts.operation.buffer.BufferOp;
import org.locationtech.jts.operation.buffer.BufferParameters;
import org.locationtech.jts.operation.union.UnaryUnionOp;
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
        private final double MIN_MU_THRESHOLD = 0.76;

        // WGS84坐标系的EPSG代码，用于定义地理坐标系统
        private final String WGS84 = "EPSG:4326";

        // 地球半径（米），用于Haversine公式计算两点间距离
        private final double R = 6371000;

        // CRS缓存，避免重复解析WKT
        private final ConcurrentHashMap<String, CoordinateReferenceSystem> crsCache = new ConcurrentHashMap<>();
        // 变换缓存，避免重复构建同一分带的投影转换
        private final ConcurrentHashMap<String, MathTransform> txCache = new ConcurrentHashMap<>();

        // GeometryFactory缓存，避免重复创建
        private final GeometryFactory geometryFactory = new GeometryFactory();

        // 默认轮廓返回的最多多边形数量（TopN）
        private final int DEFAULT_MAX_OUTLINE_SEGMENTS = 10;

        private final int DEFAULT_UNION_GROUP_SIZE = 600; // 分组合并组大小：600，避免 GeometryCollection 一次性过大
        private final int DEFAULT_BUFFER_QUADRANT = 4; // 圆近似细分：4，提升性能，误差满足道路宽度场景
        private final int DEFAULT_BUCKET_TARGET = 300; // 目标桶数：300，平衡桶内复杂度与最终合并规模
        private final int DEFAULT_BUCKET_CELL_MIN_FACTOR = 5; // cellSize 下限系数：5，避免过多小桶导致合并碎片化
        private final int DEFAULT_BUCKET_CELL_MAX_FACTOR = 36; // cellSize 上限系数：36

        private final double MIN_COMPACTNESS = 0.12; // 紧致度阈值，过滤细长道路型轮廓
        private final double MAX_ASPECT_RATIO = 8.0; // 长宽比阈值，识别道路形（越大越细长）
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

    // ===== 高斯-克吕格（6°分带）投影辅助方法：按首点经度自动分带 =====
    private int gkZoneFromLon(double lon) {
        // 6°分带：1区中心经线为3°，中心经线 L0 = 6*zone - 3
        int zone = (int) Math.floor(lon / 6.0) + 1;
        // 保护性限制（中国范围）：13–23区
        if (zone < 1)
            zone = 1;
        return zone;
    }

    private double gkCentralMeridian(int zone) {
        return zone * 6.0 - 3.0;
    }

    private double gkFalseEasting(int zone) {
        // 通用约定：假东距 = 区号 * 1,000,000 + 500,000（不会影响形状，仅影响坐标值）
        return zone * 1_000_000.0 + 500_000.0;
    }

    // 依据东距反推分带（近似），用于从高斯坐标推断区号
    private int gkZoneFromEasting(double easting) {
        int zone = (int) Math.round((easting - 500_000.0) / 1_000_000.0);
        if (zone < 1)
            zone = 1;
        if (zone > 60)
            zone = 60;
        return zone;
    }

    // 从几何的质心东距推断分带
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
    public List<TrackPoint> filterAndSortTrackPoints(List<TrackPoint> seg) {
        if (seg == null)
            return java.util.Collections.emptyList();
        java.util.Comparator<java.time.LocalDateTime> cmp = java.util.Comparator
                .nullsLast(java.util.Comparator.naturalOrder());
        return seg.stream()
                .filter(p -> p.getTime() != null)
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
        List<TrackPoint> sortedSeg = filterAndSortTrackPoints(seg);
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
            Geometry result = buildOutlineBySimpleBuffers(sortedSeg, widthM);
            long buildTime = System.currentTimeMillis() - buildStartTime;

            // 基于首点经度构建高斯-克吕格米制投影（用于面积与形态计算）
            double originLonArea = sortedSeg.get(0).getLon();
            MathTransform txWgsToGkArea = getTxWgsToGaussByLon(originLonArea);

            // 基于宽度的最小面积阈值（直接使用宽度，不再依赖配置因子）
            double minAreaThresholdM2 = Math.PI * widthM * widthM;
            if (result instanceof MultiPolygon) {
                MultiPolygon mp = (MultiPolygon) result;
                List<Polygon> keep = new ArrayList<>();
                for (int i = 0; i < mp.getNumGeometries(); i++) {
                    Polygon poly = (Polygon) mp.getGeometryN(i);
                    Geometry projPoly = JTS.transform(poly, txWgsToGkArea);
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

        // 基于首点经度构建高斯-克吕格投影转换（6°分带）
        double originLon = seg.get(0).getLon();
        MathTransform txWgsToGk = getTxWgsToGaussByLon(originLon);

        for (int i = 0; i < seg.size(); i++) {
            long startPoint = System.currentTimeMillis();

            TrackPoint point = seg.get(i);
            Coordinate coord = new Coordinate(point.getLon(), point.getLat());
            Geometry pointGeom = config.geometryFactory.createPoint(coord);

            // 转换到高斯-克吕格投影坐标系（米制）
            long startConvert = System.currentTimeMillis();
            Geometry projPoint = JTS.transform(pointGeom, txWgsToGk);
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

        // 分桶：在高斯-克吕格投影平面按网格将缓冲分组，降低相互参与的多边形数量
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

        // 转换回WGS84坐标系（使用高斯-克吕格反转换）
        long startBackConvert = System.currentTimeMillis();
        MathTransform txGkToWgs = getTxGaussToWgsByLon(originLon);
        Geometry result = JTS.transform(union, txGkToWgs);
        long backConvertTime = System.currentTimeMillis() - startBackConvert;

        long endTime = System.currentTimeMillis();
        log.trace("[buildOutlineBySimpleBuffers] 处理完成，总计耗时: {}ms (转换:{}ms, 缓冲:{}ms, 合并:{}ms, 回转:{}ms)",
                endTime - startTime, convertTime, bufferTime, unionTime, backConvertTime);

        return result;
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

    // 形状过滤：按紧致度+长宽比移除细长道路型多边形；当全部被筛掉时保留紧致度最高者
    private Geometry removeElongatedPolygons(Geometry geometry, double minCompactness, double maxAspectRatio) {
        if (geometry == null)
            return null;
        GeometryFactory gf = geometry.getFactory();
        if (geometry instanceof Polygon) {
            Polygon p = (Polygon) geometry;
            double comp = compactness(p);
            double ar = aspectRatioMeters(p);
            if (Double.isFinite(comp) && Double.isFinite(ar) && comp < minCompactness && ar > maxAspectRatio) {
                return gf.createMultiPolygon(new Polygon[0]);
            }
            return p;
        } else if (geometry instanceof MultiPolygon) {
            MultiPolygon mp = (MultiPolygon) geometry;
            java.util.List<Polygon> kept = new java.util.ArrayList<>();
            Polygon best = null;
            double bestComp = -1.0;
            for (int i = 0; i < mp.getNumGeometries(); i++) {
                Polygon p = (Polygon) mp.getGeometryN(i);
                double comp = compactness(p);
                double ar = aspectRatioMeters(p);
                if (comp >= minCompactness || ar <= maxAspectRatio) {
                    kept.add(p);
                }
                if (comp > bestComp) {
                    bestComp = comp;
                    best = p;
                }
            }
            if (kept.isEmpty()) {
                return (best != null) ? best : gf.createMultiPolygon(new Polygon[0]);
            }
            if (kept.size() == 1) {
                return kept.get(0);
            }
            return gf.createMultiPolygon(kept.toArray(new Polygon[0]));
        }
        return geometry;
    }

    // 紧致度：4πA / P²（A为球面面积平方米，P为周长米数）
    private double compactness(Polygon polygon) {
        if (polygon == null)
            return 0.0;
        double area = Math.abs(ringArea(polygon.getExteriorRing()));
        for (int i = 0; i < polygon.getNumInteriorRing(); i++) {
            area -= Math.abs(ringArea(polygon.getInteriorRingN(i)));
        }
        if (area <= 0)
            return 0.0;
        double perimeter = ringPerimeterMeters(polygon.getExteriorRing());
        for (int i = 0; i < polygon.getNumInteriorRing(); i++) {
            perimeter += ringPerimeterMeters(polygon.getInteriorRingN(i));
        }
        if (perimeter <= 0)
            return 0.0;
        return 4.0 * Math.PI * area / (perimeter * perimeter);
    }

    // 环周长（米），按haversine累加
    private double ringPerimeterMeters(org.locationtech.jts.geom.LineString ring) {
        org.locationtech.jts.geom.Coordinate[] coords = ring.getCoordinates();
        if (coords == null || coords.length < 2)
            return 0.0;
        double sum = 0.0;
        for (int i = 1; i < coords.length; i++) {
            CoordinatePoint a = new CoordinatePoint(coords[i - 1].x, coords[i - 1].y);
            CoordinatePoint b = new CoordinatePoint(coords[i].x, coords[i].y);
            sum += haversine(a, b);
        }
        return sum;
    }

    // 轴对齐包围盒的长宽比（米），用于识别道路型细长形状
    private double aspectRatioMeters(Polygon polygon) {
        Envelope env = polygon.getEnvelopeInternal();
        double midLat = (env.getMinY() + env.getMaxY()) / 2.0;
        double midLon = (env.getMinX() + env.getMaxX()) / 2.0;
        double widthM = haversine(new CoordinatePoint(env.getMinX(), midLat),
                new CoordinatePoint(env.getMaxX(), midLat));
        double heightM = haversine(new CoordinatePoint(midLon, env.getMinY()),
                new CoordinatePoint(midLon, env.getMaxY()));
        double a = Math.max(widthM, heightM);
        double b = Math.max(1e-6, Math.min(widthM, heightM));
        return a / b;
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
    /**
     * 计算两点间的大圆距离（Haversine公式）
     *
     * @return 距离（米），两点间的地理距离
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
            if (wgs instanceof org.locationtech.jts.geom.Polygon) {
                org.locationtech.jts.geom.Polygon p = (org.locationtech.jts.geom.Polygon) wgs;
                double areaOuter = Math.abs(ringArea(p.getExteriorRing()));
                double holesArea = 0.0;
                for (int i = 0; i < p.getNumInteriorRing(); i++) {
                    holesArea += Math.abs(ringArea(p.getInteriorRingN(i)));
                }
                areaSqm = areaOuter - holesArea;
            } else if (wgs instanceof org.locationtech.jts.geom.MultiPolygon) {
                org.locationtech.jts.geom.MultiPolygon mp = (org.locationtech.jts.geom.MultiPolygon) wgs;
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
        if (wkt == null || wkt.trim().isEmpty()) {
            return 0.0;
        }
        try {
            WKTReader reader = new WKTReader(config.geometryFactory);
            Geometry wgs = reader.read(wkt);
            double areaSqm = 0.0;
            if (wgs instanceof org.locationtech.jts.geom.Polygon) {
                org.locationtech.jts.geom.Polygon p = (org.locationtech.jts.geom.Polygon) wgs;
                double areaOuter = Math.abs(ringArea(p.getExteriorRing()));
                double holesArea = 0.0;
                for (int i = 0; i < p.getNumInteriorRing(); i++) {
                    holesArea += Math.abs(ringArea(p.getInteriorRingN(i)));
                }
                areaSqm = areaOuter - holesArea;
            } else if (wgs instanceof org.locationtech.jts.geom.MultiPolygon) {
                org.locationtech.jts.geom.MultiPolygon mp = (org.locationtech.jts.geom.MultiPolygon) wgs;
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
            return Math.round((areaSqm / 666.6667) * 10000.0) / 10000.0;
        } catch (org.locationtech.jts.io.ParseException e) {
            throw new RuntimeException("解析WKT字符串时出错: " + e.getMessage(), e);
        }
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
        } catch (org.locationtech.jts.io.ParseException e) {
            throw new RuntimeException("解析WKT字符串时出错: " + e.getMessage(), e);
        } catch (Exception e) {
            throw new RuntimeException("WGS84→Gauss投影转换失败: " + e.getMessage(), e);
        }
    }

    /**
     * 计算两个WKT（WGS84）的相交部分，返回相交WKT与亩数。
     * 如果无相交或相交结果为空几何，则返回 wkt=null，mu=0。
     * 支持 POLYGON 与 MULTIPOLYGON 输入，输出为相交几何的 WKT 字符串。
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
     * 基于轨迹点与总宽度生成轮廓，不会进行拆分操作
     * 
     * @param seg 轨迹点列表
     * @param totalWidthM 总宽度（米）
     * @return 轮廓部分
     * @throws Exception 异常
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
        java.time.LocalDateTime startTime = null;
        java.time.LocalDateTime endTime = null;
        for (TrackPoint p : points) {
            if (p.getTime() != null) {
                if (startTime == null)
                    startTime = p.getTime();
                endTime = p.getTime();
            }
        }

        // 使用点圆缓冲 + 方格分桶算法获取轮廓（WGS84坐标系）
        long tBuildStart = System.currentTimeMillis();
        Geometry outline = buildOutlineBySimpleBuffers(points, widthM);
        long tBuildEnd = System.currentTimeMillis();
        log.trace("[getOutline] 轮廓构建完成，耗时: {}ms", (tBuildEnd - tBuildStart));

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

                double r = Math.max(1.0, widthM); // 扩张半径：用单侧宽度作为桥接尺度
                Geometry dilated = BufferOp.bufferOp(proj, r, bp);
                Geometry hull = dilated.convexHull();
                double shrink = Math.max(0.5, widthM * 0.9); // 收缩半径：略小于扩张，保留凹特征
                Geometry eroded = BufferOp.bufferOp(hull, -shrink, bp);
                Geometry concaveWgs = JTS.transform(eroded, txGkToWgsOutline);

                if (concaveWgs instanceof Polygon) {
                    outline = (Polygon) concaveWgs;
                } else if (concaveWgs instanceof MultiPolygon && concaveWgs.getNumGeometries() == 1) {
                    outline = (Polygon) ((MultiPolygon) concaveWgs).getGeometryN(0);
                } else {
                    // 兜底：再做一次轻微闭运算后选单个外轮廓
                    Geometry closed = BufferOp.bufferOp(concaveWgs, Math.max(0.5, widthM * 0.25), bp);
                    closed = BufferOp.bufferOp(closed, -Math.max(0.5, widthM * 0.25), bp);
                    if (closed instanceof Polygon) {
                        outline = (Polygon) closed;
                    } else if (closed instanceof MultiPolygon) {
                        // 若仍为多面，合为外包凹形（最后兜底保证 Polygon）
                        Geometry finalHull = UnaryUnionOp.union(closed).convexHull();
                        if (finalHull instanceof Polygon) {
                            outline = (Polygon) finalHull;
                        } else if (finalHull instanceof MultiPolygon && finalHull.getNumGeometries() == 1) {
                            outline = (Polygon) ((MultiPolygon) finalHull).getGeometryN(0);
                        } else {
                            // 极端兜底：矩形外包
                            Geometry env = closed.getEnvelope();
                            outline = (env instanceof Polygon)
                                    ? (Polygon) env
                                    : (Polygon) ((MultiPolygon) env).getGeometryN(0);
                        }
                    }
                }
                log.trace("[getOutline] 已将MultiPolygon转换为单个Polygon（凹包近似）");
            } catch (Exception ex) {
                log.warn("[getOutline] 凹包转换失败，使用凸包兜底: {}", ex.getMessage());
                Geometry hull = outline.convexHull();
                outline = (hull instanceof Polygon)
                        ? (Polygon) hull
                        : (Polygon) ((MultiPolygon) hull).getGeometryN(0);
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
     * @return            轮廓几何，可能为 `Polygon` 或 `MultiPolygon`（裁剪后最多保留前 N 个）
     * @throws Exception  坐标转换或几何运算异常
     */
    public SplitRoadResult splitRoad(List<TrackPoint> seg, double totalWidthM) throws Exception {
        return splitRoad(seg, totalWidthM, config.DEFAULT_MAX_OUTLINE_SEGMENTS);
    }

    /**
     * 基于轨迹点与总宽度生成轮廓，并按面积保留不超过 maxSegments 的最大区块。
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
        int partsAfterAreaFilter = (trimmed instanceof MultiPolygon) ? trimmed.getNumGeometries() : 1;
        int removedByArea = partsAfterKeep - partsAfterAreaFilter;
        // 形状过滤：紧致度+长宽比，避免道路型细长形状
        trimmed = removeElongatedPolygons(trimmed, config.MIN_COMPACTNESS, config.MAX_ASPECT_RATIO);
        int partsAfterFilter = (trimmed instanceof MultiPolygon) ? trimmed.getNumGeometries() : 1;
        int removedByShape = partsAfterAreaFilter - partsAfterFilter;
        long tTrimEnd = System.currentTimeMillis();
        log.trace("[splitRoad] 裁剪+面积+形状过滤完成 type={} parts={} 移除面积={} 移除形状={} 阈值紧致度={} 阈值长宽比={} 耗时={}ms",
                trimmed.getGeometryType(), partsAfterFilter, removedByArea, removedByShape,
                config.MIN_COMPACTNESS, config.MAX_ASPECT_RATIO, (tTrimEnd - tTrimStart));

        // 最终返回日志：当少于上限时说明原因，否则保持原样
        if (partsAfterFilter < limit) {
            String reason;
            if (removedByShape > 0) {
                reason = "形状过滤(紧致度/长宽比)";
            } else if (partsAfterFilter != partsAfterKeep) {
                reason = "小面积过滤";
            } else if (partsOutline < limit) {
                reason = "原始区块数小于上限";
            } else {
                reason = "面积裁剪后区块数不超过上限";
            }
            log.debug("[splitRoad] 返回结果 type={} parts={} (上限={}，原因：{}；小面积移除={}；形状移除={})", trimmed.getGeometryType(),
                    partsAfterFilter, limit, reason, removedByArea, removedByShape);
        } else {
            log.debug("[splitRoad] 返回结果 type={} parts={}", trimmed.getGeometryType(), partsAfterFilter);
        }

        // 准备有效的轨迹点（与轮廓构建一致的过滤逻辑）
        long tFilterStart = System.currentTimeMillis();
        List<TrackPoint> sortedSeg = filterAndSortTrackPoints(seg);
        long tFilterEnd = System.currentTimeMillis();
        log.trace("[splitRoad] 点过滤+排序完成 有效点={} 耗时={}ms", sortedSeg.size(), (tFilterEnd - tFilterStart));

        java.util.List<OutlinePart> parts = new java.util.ArrayList<>();

        if (trimmed instanceof Polygon) {
            log.trace("[splitRoad] 进入Polygon分支");
            Polygon poly = (Polygon) trimmed;
            PreparedGeometry preparedPoly = PreparedGeometryFactory.prepare(poly);
            Envelope env = poly.getEnvelopeInternal();
            // 找出属于该区块的点（点在多边形内），并统计“会话段”（连续在内的首段）
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
            log.trace("[splitRoad] Polygon 会话段统计完成 sessions={} 耗时={}ms", sessions.size(), (tWithinEnd - tWithinStart));

            // 展平所有在内点集合（几何与点集保持原样），但时间口径采用首段
            java.util.List<TrackPoint> inPoly = sessions.stream()
                    .flatMap(java.util.List::stream)
                    .collect(java.util.stream.Collectors.toList());

            java.time.LocalDateTime start = sortedSeg.isEmpty() ? null : sortedSeg.get(0).getTime();
            java.time.LocalDateTime end = sortedSeg.isEmpty() ? null : sortedSeg.get(sortedSeg.size() - 1).getTime();
            log.trace("[splitRoad] Polygon 跨度(插值) 时间范围 start={} end={} inPoly={} sessions={}", start, end,
                    inPoly.size(), sessions.size());

            long tMuStart = System.currentTimeMillis();
            double mu = calcMu(poly);
            long tMuEnd = System.currentTimeMillis();
            log.trace("[splitRoad] Polygon mu计算完成 值={} 耗时={}ms", mu, (tMuEnd - tMuStart));

            long tWktStart = System.currentTimeMillis();
            String wkt = toWkt(poly);
            long tWktEnd = System.currentTimeMillis();
            log.trace("[splitRoad] Polygon WKT生成完成 长度={} 耗时={}ms", wkt.length(), (tWktEnd - tWktStart));

            // 传入轮廓内轨迹点（已按时间升序），便于后续分析
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

                java.time.LocalDateTime start = sortedSeg.isEmpty() ? null : sortedSeg.get(0).getTime();
                java.time.LocalDateTime end = sortedSeg.isEmpty() ? null
                        : sortedSeg.get(sortedSeg.size() - 1).getTime();
                log.trace("[splitRoad] Part#{} 跨度(简化) 时间范围 start={} end={} (sessions={})", i, start, end,
                        sessions.size());

                long tMuStart = System.currentTimeMillis();
                double mu = calcMu(poly);
                long tMuEnd = System.currentTimeMillis();
                log.trace("[splitRoad] Part#{} mu计算完成 值={} 耗时={}ms", i, mu, (tMuEnd - tMuStart));

                long tWktStart = System.currentTimeMillis();
                String wkt = toWkt(poly);
                long tWktEnd = System.currentTimeMillis();
                log.trace("[splitRoad] Part#{} WKT生成完成 长度={} 耗时={}ms", i, wkt.length(), (tWktEnd - tWktStart));

                // 展平会话段，形成该区块内的轨迹点集合
                java.util.List<TrackPoint> inPoly = sessions.stream()
                        .flatMap(java.util.List::stream)
                        .collect(java.util.stream.Collectors.toList());
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
        log.debug("[splitRoad] 总耗时={}ms", (tTotalEnd - t0));
        return new SplitRoadResult(trimmed, parts, outlineWkt).setTotalWidthM(totalWidthM);
    }

}
