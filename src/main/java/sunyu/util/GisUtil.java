package sunyu.util;

import java.time.Duration;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.concurrent.ConcurrentHashMap;
import java.util.stream.Collectors;

import org.geotools.geometry.jts.JTS;
import org.geotools.geometry.jts.JTSFactoryFinder;
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
import org.locationtech.jts.geom.prep.PreparedGeometry;
import org.locationtech.jts.geom.prep.PreparedGeometryFactory;
import org.locationtech.jts.io.ParseException;
import org.locationtech.jts.io.WKTReader;
import org.locationtech.jts.simplify.DouglasPeuckerSimplifier;
import org.opengis.referencing.crs.CoordinateReferenceSystem;
import org.opengis.referencing.operation.MathTransform;

import cn.hutool.core.collection.CollUtil;
import cn.hutool.log.Log;
import cn.hutool.log.LogFactory;
import sunyu.util.pojo.CoordinatePoint;
import sunyu.util.pojo.OutlinePart;
import sunyu.util.pojo.SplitRoadResult;
import sunyu.util.pojo.TrackPoint;
import sunyu.util.pojo.WktIntersectionResult;

/**
 * GIS工具类，提供全面的地理空间数据处理功能，包括坐标转换、几何计算、轨迹处理和面积计算。
 * <p>
 * 主要功能：
 * <ul>
 *   <li><strong>坐标转换</strong>：WGS84与高斯-克吕格投影坐标系之间的相互转换，支持自动投影带选择</li>
 *   <li><strong>几何计算</strong>：点与几何图形关系判断、两点间球面距离计算、几何图形相交分析等</li>
 *   <li><strong>轨迹处理</strong>：轨迹点抽稀、过滤、轨迹分段、道路拆分等农业轨迹处理功能</li>
 *   <li><strong>面积计算</strong>：支持球面面积计算，使用Haversine公式确保高精度，与Turf.js算法保持一致</li>
 *   <li><strong>性能优化</strong>：包含空间索引、PreparedGeometry、分块处理等多种优化策略</li>
 * </ul>
 * </p>
 * <p>
 * 使用说明：
 * <ul>
 *   <li>本类实现了AutoCloseable接口，建议使用try-with-resources语句管理资源</li>
 *   <li>采用Builder模式创建实例，可自定义配置参数如地球半径、投影精度等</li>
 *   <li>所有公开方法的输入和输出坐标均采用WGS84坐标系（经纬度）</li>
 *   <li>内部计算会根据需要转换至高斯投影以提高精度，特别适合需要精确面积计算的场景</li>
 *   <li>针对大数据量场景进行了性能优化，支持大量轨迹点的高效处理</li>
 * </ul>
 * </p>
 * <p>
 * 适用场景：
 * <ul>
 *   <li>农业轨迹数据分析与处理</li>
 *   <li>地理围栏与空间关系判断</li>
 *   <li>地理信息系统数据处理</li>
 *   <li>高精度面积计算与测量</li>
 *   <li>轨迹点抽稀与优化</li>
 * </ul>
 * </p>
 *
 * @author SunYu
 * @since 1.0
 * @see AutoCloseable
 */
public class GisUtil implements AutoCloseable {
    private final Log log = LogFactory.get();
    // 配置参数，包含各种常量和默认值
    private final Config config;

    /**
     * 创建GisUtil构建器实例，采用Builder模式进行对象创建。
     * <p>
     * 这是获取GisUtil实例的唯一入口，确保所有实例都经过正确配置和初始化。
     * Builder模式提供了灵活的参数配置选项，同时保证对象的一致性和完整性。
     * </p>
     * 
     * @return GisUtil构建器实例，用于配置和创建GisUtil对象
     * @see Builder
     */
    public static Builder builder() {
        return new Builder();
    }

    /**
     * 私有构造函数，通过Builder模式创建GisUtil实例。
     * <p>
     * 确保所有GisUtil实例只能通过builder()方法获取，保证配置的一致性和完整性。
     * 构造函数内部初始化必要的日志记录和配置引用。
     * </p>
     * 
     * @param config 内部配置对象，包含GIS处理所需的常量、默认值和缓存实例
     */
    private GisUtil(Config config) {
        log.info("[构建{}] 开始", this.getClass().getSimpleName());
        this.config = config;
        log.info("[构建{}] 结束", this.getClass().getSimpleName());
    }

    /**
     * 内部配置类，封装GisUtil所需的所有配置参数和共享资源。
     * <p>
     * 该类负责管理GIS处理过程中的核心配置，包括几何工厂、坐标参考系统、
     * 线程安全的缓存机制以及各种常量定义。所有配置参数都设计为不可变，
     * 确保在多线程环境下的安全访问。
     * </p>
     */
    private static class Config {
        /** 几何工厂，用于创建各种几何对象 */
        private final GeometryFactory geometryFactory = JTSFactoryFinder.getGeometryFactory();

        /** 空几何集合，用于表示无效或空的几何结果 */
        private final Geometry EMPTYGEOM = geometryFactory.createGeometryCollection();

        /** WGS84坐标参考系统，使用默认地理CRS，避免重复创建 */
        private final CoordinateReferenceSystem WGS84_CRS = DefaultGeographicCRS.WGS84;

        /**
         * 高斯投影CRS缓存（线程安全）
         * <p>
         * 键格式："zone_falseEasting_centralMeridian"
         * 用于缓存不同投影参数的坐标参考系统，避免重复创建提高性能
         * </p>
         */
        private final ConcurrentHashMap<String, CoordinateReferenceSystem> gaussCRSCache = new ConcurrentHashMap<>();

        /**
         * WGS84到高斯投影的坐标转换缓存（线程安全）
         * <p>
         * 键格式："zone_falseEasting_centralMeridian"
         * 用于缓存不同参数的坐标转换对象，显著提升重复坐标转换的性能
         * </p>
         */
        private final ConcurrentHashMap<String, MathTransform> wgs84ToGaussTransformCache = new ConcurrentHashMap<>();

        /**
         * 高斯投影到WGS84的坐标转换缓存（线程安全）
         * <p>
         * 键格式："zone_falseEasting_centralMeridian"
         * 用于缓存不同参数的坐标转换对象，显著提升重复坐标转换的性能
         * </p>
         */
        private final ConcurrentHashMap<String, MathTransform> gaussToWgs84TransformCache = new ConcurrentHashMap<>();

        /** WGS84椭球长半轴（米），与Turf.js保持一致，用于计算球面距离 */
        private final double EARTH_RADIUS = 6378137.0;

        /** 最小作业幅宽阈值（米），低于此值认为参数无效，用于农业轨迹处理 */
        private final double MIN_WORKING_WIDTH_M = 1.0;

        /** 最大返回多边形数量 */
        private final int MAX_GEOMETRY = 10;

        /** 最大作业速度(km/h) */
        private final double MAX_WORK_SPEED_1s = 8;
        private final double MAX_WORK_SPEED_10s = 13;

        /** 最小作业点阈值 */
        private final int MIN_WORK_POINTS_1s = 60;
        private final int MIN_WORK_POINTS_10s = 30;

        /** 最小亩数阈值 */
        private final double MIN_MU_1s = 0.3;
        private final double MIN_MU_10s = 0.5;
    }

    /**
     * GisUtil构建器类，提供流式API用于配置和创建GisUtil实例。
     * <p>
     * 该类实现了Builder设计模式，允许在创建GisUtil对象前灵活配置各项参数。
     * 构建器模式确保对象创建的一致性和完整性，避免了构造函数参数膨胀问题。
     * </p>
     * <p>
     * 使用示例：
     * <pre>
     * try (GisUtil gisUtil = GisUtil.builder().build()) {
     *     // 使用gisUtil进行GIS操作
     * }
     * </pre>
     * </p>
     */
    public static class Builder {
        /** 配置对象，包含GisUtil所需的常量、默认值和缓存实例 */
        private Config config = new Config();

        /**
         * 构建并返回完全初始化的GisUtil实例。
         * <p>
         * 调用此方法会使用当前配置参数创建GisUtil对象，并完成所有必要的初始化。
         * 返回的对象已准备就绪，可以直接用于各种GIS操作。
         * </p>
         * 
         * @return 初始化完成的GisUtil实例
         * @see GisUtil#builder()
         */
        public GisUtil build() {
            return new GisUtil(config);
        }
    }

    /**
     * 关闭GisUtil实例，释放相关资源。
     * <p>
     * 实现AutoCloseable接口的方法，用于释放GisUtil使用的各种资源。
     * 建议使用try-with-resources语句自动调用此方法，确保资源正确释放。
     * </p>
     * 
     * @see AutoCloseable#close()
     */
    @Override
    public void close() {
        log.info("[销毁{}] 开始", this.getClass().getSimpleName());
        // 可以在此处添加资源清理逻辑，如清空缓存等
        log.info("[销毁{}] 结束", this.getClass().getSimpleName());
    }

    /**
     * 使用Douglas-Peucker算法对轨迹点进行抽稀
     * 
     * @param points    原始轨迹点列表
     * @param tolerance 容差（米）
     * @return 抽稀后的轨迹点列表
     */
    private List<TrackPoint> simplifyTrackPoints(List<TrackPoint> points, double tolerance) {
        if (points.size() <= 2) {
            return new ArrayList<>(points);
        }

        // 转换为Coordinate数组用于Douglas-Peucker算法
        Coordinate[] coords = new Coordinate[points.size()];
        for (int i = 0; i < points.size(); i++) {
            TrackPoint p = points.get(i);
            coords[i] = new Coordinate(p.getLon(), p.getLat());
        }

        // 执行Douglas-Peucker抽稀
        // 先创建LineString几何对象，再应用简化算法
        GeometryFactory factory = new GeometryFactory();
        LineString line = factory.createLineString(coords);
        Geometry simplifiedGeometry = DouglasPeuckerSimplifier.simplify(line, tolerance);
        Coordinate[] simplifiedCoords = simplifiedGeometry.getCoordinates();

        // 将结果映射回TrackPoint列表
        List<TrackPoint> result = new ArrayList<>(simplifiedCoords.length);
        Map<Coordinate, TrackPoint> coordToPointMap = new HashMap<>();

        // 构建坐标到点的映射，优先保留时间信息
        for (TrackPoint point : points) {
            Coordinate coord = new Coordinate(point.getLon(), point.getLat());
            coordToPointMap.put(coord, point);
        }

        // 匹配抽稀后的坐标对应的原始点
        for (Coordinate coord : simplifiedCoords) {
            // 精确匹配
            TrackPoint exactMatch = coordToPointMap.get(coord);
            if (exactMatch != null) {
                result.add(exactMatch);
            } else {
                // 如果没有精确匹配，查找最接近的点
                TrackPoint closest = findClosestPoint(coord, points);
                if (closest != null) {
                    result.add(closest);
                }
            }
        }

        return result;
    }

    /**
     * 查找最接近给定坐标的轨迹点
     */
    private TrackPoint findClosestPoint(Coordinate coord, List<TrackPoint> points) {
        TrackPoint closest = null;
        double minDist = Double.MAX_VALUE;

        for (TrackPoint point : points) {
            // 计算欧几里得距离
            double dx = coord.x - point.getLon();
            double dy = coord.y - point.getLat();
            double dist = Math.sqrt(dx * dx + dy * dy);
            if (dist < minDist) {
                minDist = dist;
                closest = point;
            }
        }

        return closest;
    }

    /**
     * 分块处理大型轨迹段 - 高级优化版本
     * 关键优化：
     * 1. 动态调整块大小和重叠区域，进一步减少几何图形复杂度
     * 2. 优化并行配置和分块策略，减少线程创建开销
     * 3. 提升buffer计算效率，添加线段简化预处理
     * 4. 优化内存使用，减少中间对象创建
     * 
     * @param points      轨迹点列表
     * @param bufferWidth 缓冲区宽度（米）
     * @return 合并后的几何图形
     */
    private Geometry processLargeSegmentInChunks(List<TrackPoint> points, double bufferWidth) {
        long startTime = System.currentTimeMillis();
        log.debug("开始分块处理大型轨迹段，点数: {}", points.size());

        // 简化分块逻辑：使用固定的较大块大小，减少分块数
        int totalPoints = points.size();
        int chunkSize = Math.max(1000, Math.min(totalPoints, 2000)); // 使用更大的块大小减少分块数
        int overlapSize = 50; // 使用固定的重叠区域确保连续性

        // 计算块数并生成块起始索引
        int chunksCount = (int) Math.ceil((double) totalPoints / (chunkSize - overlapSize));
        List<Integer> chunkStartIndices = new ArrayList<>(chunksCount);

        for (int i = 0; i < chunksCount; i++) {
            int startIndex = i * (chunkSize - overlapSize);
            // 确保最后一个块的完整性
            if (i == chunksCount - 1 && startIndex + chunkSize > totalPoints) {
                startIndex = Math.max(0, totalPoints - chunkSize);
            }
            chunkStartIndices.add(startIndex);
        }

        log.debug("共分成 {} 个块进行处理，块大小: {}, 重叠区域: {}",
                chunkStartIndices.size(), chunkSize, overlapSize);

        // 直接使用顺序处理所有块
        List<Geometry> chunkGeometries = new ArrayList<>();
        for (Integer startIndex : chunkStartIndices) {
            Geometry geom = processChunk(points, startIndex, chunkSize, totalPoints, bufferWidth);
            if (geom != null && !geom.isEmpty()) {
                chunkGeometries.add(geom);
            }
        }

        log.debug("所有块处理完成，准备合并几何图形，耗时: {}ms", System.currentTimeMillis() - startTime);

        // 直接合并几何图形，移除中间预处理步骤
        Geometry mergedGeometry = mergeGeometriesRecursively(chunkGeometries);

        // 简化后处理：只进行基本的有效性检查和修复，避免过度处理
        if (!mergedGeometry.isValid()) {
            mergedGeometry = mergedGeometry.buffer(0); // 修复无效几何图形
        }

        log.debug("分块处理完成，合并后几何类型: {}, 总耗时: {}ms",
                mergedGeometry.getGeometryType(), System.currentTimeMillis() - startTime);
        return mergedGeometry;
    }

    /**
     * 处理单个轨迹块，提取为单独方法以提高可读性和可维护性
     */
    private Geometry processChunk(List<TrackPoint> points, int startIndex, int chunkSize,
            int totalPoints, double bufferWidth) {
        int end = Math.min(startIndex + chunkSize, totalPoints);
        List<TrackPoint> chunk = points.subList(startIndex, end);

        long chunkStartTime = System.currentTimeMillis();

        // 优化：直接使用数组而不是中间集合
        int chunkLength = chunk.size();
        Coordinate[] coords = new Coordinate[chunkLength];
        for (int j = 0; j < chunkLength; j++) {
            TrackPoint p = chunk.get(j);
            coords[j] = new Coordinate(p.getLon(), p.getLat());
        }

        LineString chunkLine = config.geometryFactory.createLineString(coords);

        // 优化4: 进一步优化buffer参数，使用最小必要的精度
        int quadrantSegments = 4; // 更低的值，显著提升性能，在保证准确性的前提下
        int endCapStyle = 1; // 1=LINEAREND，高效的样式

        // 计算buffer
        Geometry chunkBuffer = chunkLine.buffer(bufferWidth, quadrantSegments, endCapStyle);

        // 可选：简化几何图形，进一步减少复杂度（如果buffer后图形仍然复杂）
        if (chunkBuffer.getNumPoints() > 1000) {
            chunkBuffer = DouglasPeuckerSimplifier.simplify(chunkBuffer, 0.00001); // 轻微简化，保持形状
        }

        log.debug("处理块: {}-{}, 点数: {}, buffer几何类型: {}, 耗时: {}ms",
                startIndex, end, chunkLength, chunkBuffer.getGeometryType(),
                System.currentTimeMillis() - chunkStartTime);

        return chunkBuffer;
    }

    /**
     * 高效合并几何图形列表
     * 使用更优的合并策略，优化性能
     * 
     * @param geometries 几何图形列表
     * @return 合并后的几何图形
     */
    private Geometry mergeGeometriesRecursively(List<Geometry> geometries) {
        // 基本情况：列表为空或只有一个元素
        if (geometries.isEmpty()) {
            return config.EMPTYGEOM;
        }
        if (geometries.size() == 1) {
            return geometries.get(0);
        }

        // 优化合并策略：使用空间索引提高合并效率
        // 1. 过滤掉空的几何图形
        List<Geometry> nonEmptyGeometries = new ArrayList<>();
        for (Geometry geom : geometries) {
            if (geom != null && !geom.isEmpty()) {
                nonEmptyGeometries.add(geom);
            }
        }

        if (nonEmptyGeometries.isEmpty()) {
            return config.EMPTYGEOM;
        }
        if (nonEmptyGeometries.size() == 1) {
            return nonEmptyGeometries.get(0);
        }

        // 2. 使用空间索引优化合并顺序
        // 优先合并相邻的几何图形，减少计算量
        Geometry result = nonEmptyGeometries.get(0);
        List<Geometry> remainingGeometries = new ArrayList<>(nonEmptyGeometries.subList(1, nonEmptyGeometries.size()));

        while (!remainingGeometries.isEmpty()) {
            boolean merged = false;
            for (int i = 0; i < remainingGeometries.size(); i++) {
                Geometry nextGeom = remainingGeometries.get(i);

                // 快速检查边界框是否相交，避免不必要的union操作
                if (result.getEnvelopeInternal().intersects(nextGeom.getEnvelopeInternal())) {
                    try {
                        Geometry newResult = result.union(nextGeom);
                        // 只有当合并成功（面积增加）时才更新结果
                        if (!newResult.isEmpty()) {
                            result = newResult;
                            remainingGeometries.remove(i);
                            merged = true;
                            break;
                        }
                    } catch (Exception e) {
                        log.trace("合并几何图形时出错: {}", e.getMessage());
                        // 合并失败时继续下一个几何图形
                    }
                }
            }

            // 如果没有找到可合并的几何图形，强制合并第一个
            if (!merged && !remainingGeometries.isEmpty()) {
                result = result.union(remainingGeometries.remove(0));
            }
        }

        return result;
    }

    /**
     * 计算线环的球面面积（平方米）
     * <p>
     * 该方法使用球面多边形面积计算公式，精确考虑地球曲率，适用于任意区域大小的面积计算。
     * 算法与Turf.js保持一致，确保跨平台结果兼容性。公式为：A = R² × |Σ(λi+1 - λi) × sin((φi+1 + φi)/2)|
     * 其中R为地球半径，λ为经度（弧度），φ为纬度（弧度）。
     * </p>
     * 
     * @param wgs84Ring WGS84坐标系下的线环（单位：度），必须是闭合的环
     * @return 球面面积（平方米），计算失败返回0.0
     */
    private double calculateRingSphericalArea(LineString wgs84Ring) {
        // 参数有效性检查
        if (wgs84Ring == null || wgs84Ring.isEmpty()) {
            return 0.0;
        }

        // 获取环的所有坐标点
        org.locationtech.jts.geom.Coordinate[] coords = wgs84Ring.getCoordinates();
        // 至少需要3个点才能形成多边形
        if (coords.length < 3) {
            return 0.0;
        }

        double area = 0.0;

        // 步骤1：遍历环上的连续点对，应用球面面积公式
        for (int i = 0; i < coords.length - 1; i++) {
            // 步骤2：将经纬度转换为弧度
            double lon1 = Math.toRadians(coords[i].x);
            double lat1 = Math.toRadians(coords[i].y);
            double lon2 = Math.toRadians(coords[i + 1].x);
            double lat2 = Math.toRadians(coords[i + 1].y);

            // 步骤3：应用球面面积公式的累加项
            // 计算：(λi+1 - λi) × sin((φi+1 + φi)/2)
            area += (lon2 - lon1) * Math.sin((lat1 + lat2) / 2.0);
        }

        // 步骤4：计算最终面积 = |面积累加值| × 地球半径²
        area = Math.abs(area) * config.EARTH_RADIUS * config.EARTH_RADIUS;
        return area;
    }

    /**
     * 计算单个多边形的球面面积（平方米）
     * <p>
     * 该方法精确计算WGS84坐标系下多边形的球面面积，同时考虑了地球曲率和多边形内部孔洞。
     * 实现方式为：计算外环面积减去所有内环（孔洞）面积的总和。支持复杂多边形形状，包括带孔洞的多边形。
     * 适用于精确面积测量、土地面积计算和地理统计分析。
     * </p>
     * 
     * @param wgs84Polygon WGS84坐标系下的多边形，可包含内部孔洞
     * @return 多边形的球面面积（平方米）
     */
    private double calculatePolygonSphericalArea(Polygon wgs84Polygon) {
        // 步骤1：计算多边形外环的面积
        double exteriorArea = calculateRingSphericalArea(wgs84Polygon.getExteriorRing());
        log.trace("外环面积: {}平方米", exteriorArea);

        // 步骤2：计算所有内环（孔洞）的面积总和
        double holesArea = 0.0;
        for (int i = 0; i < wgs84Polygon.getNumInteriorRing(); i++) {
            // 计算单个内环面积
            double holeArea = calculateRingSphericalArea(wgs84Polygon.getInteriorRingN(i));
            // 累加所有孔洞面积
            holesArea += holeArea;
            log.trace("内环{}面积: {}平方米", i, holeArea);
        }

        // 步骤3：总面积 = 外环面积 - 所有孔洞面积
        double totalArea = exteriorArea - holesArea;
        log.trace("多边形总面积: {}平方米", totalArea);
        return totalArea;
    }

    /**
     * 从缓存获取或创建高斯-克吕格投影坐标参考系统（CRS）。
     * <p>
     * 该方法实现了高斯投影CRS的高效管理，使用线程安全的缓存机制避免重复创建相同参数的CRS。
     * 采用Transverse Mercator投影（横轴墨卡托投影），适用于中纬度地区的精确测量。
     * </p>
     * 
     * @param zone 高斯投影带号（1-60，对应全球6度分带）
     * @param falseEasting 假东距（米），用于避免负值坐标
     * @param centralMeridian 中央经线（度），投影带的中心线
     * @return 高斯投影坐标参考系统，创建失败返回null
     */
    private CoordinateReferenceSystem getGaussCRS(int zone, double falseEasting, double centralMeridian) {
        // 创建缓存键，使用带号、假东距和中央经线作为唯一标识
        String cacheKey = String.format("%d_%.1f_%.1f", zone, falseEasting, centralMeridian);

        // 使用computeIfAbsent实现线程安全的缓存机制，只在缓存未命中时创建新CRS
        return config.gaussCRSCache.computeIfAbsent(cacheKey, key -> {
            try {
                log.debug("创建高斯投影CRS：投影带号={}, 假东距={}, 中央经线={}", zone, falseEasting,
                        centralMeridian);
                // 定义高斯-克吕格投影坐标系 - 使用预构建的WKT模板
                // WKT格式定义了完整的坐标系统参数，包括基准面、椭球体、投影方式和参数
                String wktTemplate = "PROJCS[\"Gauss_Kruger_ZONE_" + zone + "\"," +
                        " GEOGCS[\"GCS_WGS_1984\", DATUM[\"WGS_1984\", SPHEROID[\"WGS_84\", 6378137.0, 298.257223563]], PRIMEM[\"Greenwich\", 0.0], UNIT[\"Degree\", 0.0174532925199433]], PROJECTION[\"Transverse_Mercator\"], PARAMETER[\"False_Easting\", %.6f], PARAMETER[\"False_Northing\", 0.0], PARAMETER[\"Central_Meridian\", %.6f], PARAMETER[\"Scale_Factor\", 1.0], PARAMETER[\"Latitude_Of_Origin\", 0.0], UNIT[\"Meter\", 1.0]]";
                String gaussProjString = String.format(wktTemplate, falseEasting, centralMeridian);

                // 解析WKT字符串创建CRS对象
                return CRS.parseWKT(gaussProjString);
            } catch (Exception e) {
                log.warn("创建高斯投影CRS失败：zone={}, falseEasting={}, centralMeridian={}, 错误={}",
                        zone, falseEasting, centralMeridian, e.getMessage());
                return null;
            }
        });
    }

    /**
     * 将高斯-克吕格投影坐标系下的几何图形转换为WGS84地理坐标系下的几何图形。
     * <p>
     * 该方法实现了高精度的坐标逆转换，包括：
     * <ul>
     *   <li>根据几何中心智能推断投影带号</li>
     *   <li>坐标合理性验证</li>
     *   <li>跨越多个投影带的几何图形处理策略</li>
     *   <li>线程安全的转换缓存机制</li>
     * </ul>
     * 适用于将平面坐标数据转换回地理坐标，用于全球定位和地图显示。
     * </p>
     * 
     * @param gaussGeometry 高斯投影坐标系下的几何图形（米制）
     * @return WGS84坐标系下的几何图形（经纬度），转换失败返回空几何
     */
    public Geometry toWgs84Geometry(Geometry gaussGeometry) {
        try {
            // 获取几何图形的边界信息来反推高斯投影参数
            Envelope env = gaussGeometry.getEnvelopeInternal();
            double centerX = (env.getMinX() + env.getMaxX()) / 2.0;

            log.trace("高斯投影几何到WGS84投影几何转换开始");

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
            // 计算中央经线（与wgs84PointTransformToGaussPoint方法保持一致）
            double centralMeridian = (zone - 1) * 6 - 180 + 3;
            double falseEasting = zone * 1000000.0 + 500000.0;
            log.trace("计算投影参数：zone={}, centralMeridian={}, falseEasting={}", zone, centralMeridian, falseEasting);

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

            // 从缓存获取或创建高斯投影CRS
            CoordinateReferenceSystem gaussCRS = getGaussCRS(zone, falseEasting, centralMeridian);

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
            log.trace("高斯投影几何到WGS84投影几何转换完成");
            return wgs84Geometry;
        } catch (Exception e) {
            log.warn("高斯投影几何到WGS84投影几何转换失败：错误={}", e.getMessage());
            return config.EMPTYGEOM;
        }
    }

    /**
     * 将WGS84坐标系的WKT（Well-Known Text）字符串解析为Geometry几何图形。
     * <p>
     * 该方法将标准WKT格式的几何描述文本解析为JTS几何对象，支持所有标准几何类型，
     * 如点(POINT)、线(LINESTRING)、多边形(POLYGON)、多点(MULTIPOINT)等。
     * 解析后的几何对象可直接用于GIS空间分析操作。
     * </p>
     * 
     * @param wgs84WKT WGS84坐标系下的WKT字符串，如"POINT(116.4 39.9)"或"POLYGON((...))"
     * @return WGS84坐标系的Geometry几何图形，解析失败返回空几何
     */
    public Geometry toWgs84Geometry(String wgs84WKT) {
        if (wgs84WKT == null || wgs84WKT.trim().isEmpty()) {
            log.warn("WKT字符串为空或null");
            return config.EMPTYGEOM;
        }
        try {
            Geometry geometry = new WKTReader(config.geometryFactory).read(wgs84WKT);
            log.debug("WKT字符串解析成功：几何类型={}", geometry.getGeometryType());
            return geometry;
        } catch (ParseException e) {
            log.warn("WKT字符串解析失败：{}", e.getMessage());
            return config.EMPTYGEOM;
        } catch (Exception e) {
            log.warn("WKT字符串转换几何图形失败：{}", e.getMessage());
            return config.EMPTYGEOM;
        }
    }

    /**
     * 将WGS84地理坐标系的几何图形转换为高斯-克吕格投影坐标系。
     * <p>
     * 该方法实现了高精度的坐标正转换，包括：
     * <ul>
     *   <li>全球范围支持，自动处理1-60投影带</li>
     *   <li>根据几何中心经度智能选择最合适的投影带</li>
     *   <li>输入输出坐标合理性验证</li>
     *   <li>线程安全的转换缓存机制</li>
     * </ul>
     * 高斯投影适用于距离和面积的精确计算，特别适合中纬度地区的GIS分析。
     * </p>
     * 
     * @param wgs84Geometry WGS84坐标系下的几何图形（经纬度）
     * @return 高斯投影坐标系下的几何图形（米制），转换失败返回空几何
     */
    public Geometry toGaussGeometry(Geometry wgs84Geometry) {
        try {
            if (wgs84Geometry == null || wgs84Geometry.isEmpty()) {
                return config.EMPTYGEOM;
            }

            // 获取几何图形的边界信息来确定高斯投影参数
            Envelope env = wgs84Geometry.getEnvelopeInternal();
            double centerLon = (env.getMinX() + env.getMaxX()) / 2.0;
            double centerLat = (env.getMinY() + env.getMaxY()) / 2.0;

            log.debug("WGS84几何到高斯投影几何转换开始");
            log.trace("WGS84几何边界：MinLon={}, MaxLon={}, MinLat={}, MaxLat={}",
                    env.getMinX(), env.getMaxX(), env.getMinY(), env.getMaxY());
            log.trace("几何中心坐标：centerLon={}, centerLat={}", centerLon, centerLat);

            // 验证WGS84坐标的合理性
            if (env.getMinX() < -180 || env.getMaxX() > 180 ||
                    env.getMinY() < -90 || env.getMaxY() > 90) {
                log.warn("WGS84坐标超出合理范围：MinLon={}, MaxLon={}, MinLat={}, MaxLat={}",
                        env.getMinX(), env.getMaxX(), env.getMinY(), env.getMaxY());
                return config.EMPTYGEOM;
            }

            // 计算高斯投影带号 (6度分带)
            int zone = (int) Math.floor((centerLon + 180) / 6) + 1;
            double centralMeridian = (zone - 1) * 6 - 180 + 3;
            double falseEasting = zone * 1000000.0 + 500000.0;

            // 验证投影带号的合理性
            if (zone < 1 || zone > 60) {
                log.warn("投影带号超出合理范围：zone={}，经度={}", zone, centerLon);
                return config.EMPTYGEOM;
            }

            // 从缓存获取或创建高斯投影CRS
            CoordinateReferenceSystem gaussCRS = getGaussCRS(zone, falseEasting, centralMeridian);

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

            // 执行坐标转换（正向：WGS84 -> 高斯投影）
            log.trace("zone={}, 几何点数={}", (int) Math.floor((centerLon + 180) / 6) + 1, wgs84Geometry.getNumPoints());
            Geometry gaussGeometry = JTS.transform(wgs84Geometry, transform);

            // 验证转换后的高斯投影坐标合理性
            Envelope gaussEnv = gaussGeometry.getEnvelopeInternal();
            if (gaussEnv.getMinX() < 500000 || gaussEnv.getMaxX() > 64000000 ||
                    gaussEnv.getMinY() < -10000000 || gaussEnv.getMaxY() > 10000000) {
                log.warn("转换后的高斯投影坐标超出合理范围：MinX={}, MaxX={}, MinY={}, MaxY={}",
                        gaussEnv.getMinX(), gaussEnv.getMaxX(), gaussEnv.getMinY(), gaussEnv.getMaxY());
                return config.EMPTYGEOM;
            }
            log.debug("WGS84几何到高斯投影几何转换完成");
            return gaussGeometry;
        } catch (Exception e) {
            log.warn("WGS84几何到高斯投影几何转换失败：错误={}", e.getMessage());
            return config.EMPTYGEOM;
        }
    }

    /**
     * 将单个WGS84地理坐标系的轨迹点转换为高斯-克吕格投影坐标系的轨迹点。
     * <p>
     * 该方法对单个轨迹点执行精确的坐标转换，保持时间属性不变，同时进行：
     * <ul>
     *   <li>输入坐标有效性检查</li>
     *   <li>投影带自动计算和验证</li>
     *   <li>转换结果合理性验证</li>
     *   <li>完整的错误处理和日志记录</li>
     * </ul>
     * 适用于轨迹点的单独处理和分析。
     * </p>
     * 
     * @param wgs84Point WGS84坐标系的轨迹点（经纬度）
     * @return 高斯投影坐标系的轨迹点（米制），转换失败返回null
     */
    public TrackPoint toGaussPoint(TrackPoint wgs84Point) {
        try {
            // 获取经度和纬度
            double longitude = wgs84Point.getLon();
            double latitude = wgs84Point.getLat();

            log.trace("开始WGS84到高斯投影转换 原始坐标：经度={}, 纬度={}", longitude, latitude);

            // 计算高斯投影带号 (6度分带)
            // 全球范围: 经度-180到180，对应带号1-60
            int zone = (int) Math.floor((longitude + 180) / 6) + 1;
            log.trace("计算投影带号：zone={}", zone);

            // 验证投影带号的合理性
            if (zone < 1 || zone > 60) {
                log.warn("投影带号超出合理范围：zone={}，经度={}", zone, longitude);
                return null;
            }

            // 计算中央经线
            double centralMeridian = (zone - 1) * 6 - 180 + 3;
            log.trace("计算中央经线：centralMeridian={}", centralMeridian);

            // 全球支持模式：假东距包含带号信息，便于识别投影带
            // 格式：zone × 1000000 + 500000（如49带 = 49500000米）
            double falseEasting = zone * 1000000.0 + 500000.0;
            log.trace("计算假东距：falseEasting={}", falseEasting);

            // 从缓存获取或创建高斯投影CRS
            CoordinateReferenceSystem gaussCRS = getGaussCRS(zone, falseEasting, centralMeridian);
            if (gaussCRS == null) {
                return null;
            }

            // 构建缓存key（用于transform缓存）
            String cacheKey = String.format("%d_%.1f_%.1f", zone, falseEasting, centralMeridian);

            // 从缓存获取或创建WGS84到高斯投影的坐标转换
            MathTransform transform = config.wgs84ToGaussTransformCache.computeIfAbsent(cacheKey, key -> {
                try {
                    log.debug("创建WGS84到高斯投影坐标转换器：投影带号={}, 假东距={}, 中央经线={}", zone, falseEasting,
                            centralMeridian);
                    return CRS.findMathTransform(config.WGS84_CRS, gaussCRS, true);
                } catch (Exception e) {
                    log.warn("创建坐标转换失败：zone={}, falseEasting={}, centralMeridian={}, 错误={}",
                            zone, falseEasting, centralMeridian, e.getMessage());
                    return null;
                }
            });
            if (transform == null) {
                return null;
            }

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

            log.trace("WGS84到高斯投影转换完成");
            return result;
        } catch (Exception e) {
            log.warn("WGS84到高斯投影转换失败：经度={}, 纬度={}, 错误={}", wgs84Point.getLon(), wgs84Point.getLat(), e.getMessage());
            return null;
        }
    }

    /**
     * 将WGS84地理坐标系的轨迹点列表批量转换为高斯-克吕格投影坐标系的轨迹点列表。
     * <p>
     * 该方法对轨迹点列表进行批量处理，为每个点调用toGaussPoint方法，
     * 自动过滤转换失败的点，确保返回的列表中只包含有效转换结果。
     * 适用于轨迹数据的整体处理和分析。
     * </p>
     * 
     * @param wgs84Points WGS84坐标系的轨迹点列表（经纬度）
     * @return 高斯投影坐标系的轨迹点列表（米制），保持与输入列表相同的顺序，
     *         只包含成功转换的点，可能为空列表
     */
    public List<TrackPoint> toGaussPointList(List<TrackPoint> wgs84Points) {
        List<TrackPoint> gaussPoints = new ArrayList<>();
        try {
            // 对每个轨迹点进行坐标转换
            for (TrackPoint wgs84Point : wgs84Points) {
                TrackPoint gaussPoint = toGaussPoint(wgs84Point);
                if (gaussPoint != null) {
                    gaussPoints.add(gaussPoint);
                }
            }

            log.trace("成功转换{}个轨迹点到高斯-克吕格投影坐标系", gaussPoints.size());
        } catch (Exception e) {
            log.warn("WGS84轨迹点转换为高斯投影失败: {}", e.getMessage());
        }
        return gaussPoints;
    }

    /**
     * 使用Haversine公式计算WGS84坐标系下两点之间的球面距离。
     * <p>
     * Haversine公式适用于球面上两点之间距离的精确计算，特别适合地球表面两点间的距离估算。
     * 该方法具有以下特点：
     * <ul>
     *   <li>使用地球赤道半径作为计算基准</li>
     *   <li>考虑了经纬度的弧度转换</li>
     *   <li>数值稳定性好，适用于中短距离计算</li>
     *   <li>执行效率高，无需额外的几何库支持</li>
     * </ul>
     * 适用于轨迹分析、地理围栏、距离测量等场景。
     * </p>
     * 
     * @param wgs84Point1 第一个WGS84坐标点，包含有效经纬度
     * @param wgs84Point2 第二个WGS84坐标点，包含有效经纬度
     * @return 两点之间的球面距离（米），结果为非负数
     */
    public double haversine(CoordinatePoint wgs84Point1, CoordinatePoint wgs84Point2) {
        // 将经纬度从度数转换为弧度，三角函数计算需要弧度值
        double lon1 = Math.toRadians(wgs84Point1.getLon());
        double lat1 = Math.toRadians(wgs84Point1.getLat());
        double lon2 = Math.toRadians(wgs84Point2.getLon());
        double lat2 = Math.toRadians(wgs84Point2.getLat());

        // 计算两点间经度差和纬度差
        double dlon = lon2 - lon1;
        double dlat = lat2 - lat1;

        // Haversine公式核心计算：半正矢公式
        double a = Math.sin(dlat / 2.0) * Math.sin(dlat / 2.0) +
                Math.cos(lat1) * Math.cos(lat2) *
                        Math.sin(dlon / 2.0) * Math.sin(dlon / 2.0);

        // 计算haversine函数结果，避免数值精度问题
        double c = 2.0 * Math.atan2(Math.sqrt(a), Math.sqrt(1.0 - a));

        // 最终距离 = 地球半径 * haversine函数结果
        return config.EARTH_RADIUS * c;
    }

    /**
     * 使用Haversine公式判断一个地理点是否在指定半径的圆内。
     * <p>
     * 该方法通过计算测试点与圆心之间的球面距离，并与指定半径进行比较，
     * 来确定点是否位于圆内。采用球面距离计算确保在全球范围内的准确性，
     * 特别适合地理围栏（Geo-fencing）应用场景。
     * </p>
     * 
     * @param wgs84Point 要判断的WGS84坐标点，包含有效经纬度
     * @param wgs84CenterPoint 圆的中心点（WGS84坐标）
     * @param radius 圆的半径（米），必须为非负数
     * @return 如果点到圆心的球面距离小于等于指定半径（包含边界），则返回true；否则返回false
     */
    public boolean inCircle(CoordinatePoint wgs84Point, CoordinatePoint wgs84CenterPoint, double radius) {
        // 计算两点间的球面距离（米）
        double distance = haversine(wgs84Point, wgs84CenterPoint);
        // 如果距离小于等于半径，则点在圆内（包含边界）
        return distance <= radius;
    }

    /**
     * 判断WGS84坐标系下的点是否在指定的几何图形内部。
     * <p>
     * 该方法使用JTS拓扑套件进行精确的空间关系判断，支持各种几何类型（点、线、面等）。
     * 采用covers()方法判断，确保边界上的点也被视为在几何图形内部。
     * 适用于轨迹分析、地理围栏和空间统计等场景。
     * </p>
     * 
     * @param wgs84Point 待测试的WGS84坐标点，包含有效经纬度
     * @param wgs84Geometry 用于判断的几何图形，必须是有效的JTS Geometry对象
     * @return 如果点在几何图形内部或边界上，则返回true；否则返回false
     */
    public boolean inGeometry(CoordinatePoint wgs84Point, Geometry wgs84Geometry) {
        try {
            // 创建点几何对象
            Point point = config.geometryFactory.createPoint(
                    new Coordinate(wgs84Point.getLon(), wgs84Point.getLat()));

            // 判断点是否在几何图形内或边界上
            // covers() 方法包含边界，而 contains() 不包含边界
            return wgs84Geometry.covers(point);
        } catch (Exception e) {
            log.warn("判断点是否在几何图形内失败：点[{},{}] 错误={}",
                    wgs84Point.getLon(), wgs84Point.getLat(), e.getMessage());
            return false;
        }
    }

    /**
     * 判断WGS84坐标系下的点是否位于指定的地理矩形范围内。
     * <p>
     * 该方法使用高效的边界框检查，通过比较点的经纬度与矩形边界，快速判断点是否在矩形范围内。
     * 支持任意顺序的对角点输入，内部会自动调整为正确的边界。不考虑地球曲率，适用于较小范围区域。
     * 常用于空间索引、数据分区和快速查询优化。
     * </p>
     * 
     * @param wgs84Point 待测试的WGS84坐标点，包含有效经纬度
     * @param wgs84TopLeftPoint 矩形的左上角点（或任意对角点）
     * @param wgs84BottomRightPoint 矩形的右下角点（或任意对角点）
     * @return 如果点的经纬度在矩形范围内（包含边界），则返回true；否则返回false
     */
    public boolean inRectangle(CoordinatePoint wgs84Point, CoordinatePoint wgs84TopLeftPoint,
            CoordinatePoint wgs84BottomRightPoint) {
        try {
            double pointLon = wgs84Point.getLon();
            double pointLat = wgs84Point.getLat();
            double topLeftLon = wgs84TopLeftPoint.getLon();
            double topLeftLat = wgs84TopLeftPoint.getLat();
            double bottomRightLon = wgs84BottomRightPoint.getLon();
            double bottomRightLat = wgs84BottomRightPoint.getLat();

            // 确保左上角和右下角的经纬度关系正确
            double minLon = Math.min(topLeftLon, bottomRightLon);
            double maxLon = Math.max(topLeftLon, bottomRightLon);
            double maxLat = Math.max(topLeftLat, bottomRightLat);
            double minLat = Math.min(topLeftLat, bottomRightLat);

            // 判断点是否在矩形范围内（包含边界）
            return pointLon >= minLon && pointLon <= maxLon &&
                    pointLat >= minLat && pointLat <= maxLat;
        } catch (Exception e) {
            log.warn("判断点是否在矩形内失败：点[{},{}] 左上角[{},{}] 右下角[{},{}] 错误={}",
                    wgs84Point.getLon(), wgs84Point.getLat(),
                    wgs84TopLeftPoint.getLon(), wgs84TopLeftPoint.getLat(),
                    wgs84BottomRightPoint.getLon(), wgs84BottomRightPoint.getLat(),
                    e.getMessage());
            return false;
        }
    }

    /**
     * 计算WGS84坐标系下两个几何图形的空间交集轮廓及面积。
     * <p>
     * 该方法通过以下步骤计算两个几何图形的交集：
     * <ol>
     *   <li>解析WKT字符串为JTS几何对象</li>
     *   <li>将几何图形转换到高斯投影坐标系进行高精度计算</li>
     *   <li>执行几何交集操作</li>
     *   <li>将结果转换回WGS84坐标系</li>
     *   <li>计算相交区域面积并转换为亩数</li>
     * </ol>
     * 使用高斯投影进行中间计算可显著提高面积计算精度，适用于空间分析、数据裁剪和覆盖分析。
     * </p>
     * 
     * @param wgs84WKT1 第一个WGS84几何图形的WKT字符串（如多边形、线、点等）
     * @param wgs84WKT2 第二个WGS84几何图形的WKT字符串（如多边形、线、点等）
     * @return 包含相交轮廓WKT字符串和面积（亩）的结果对象，无相交则返回空几何和零面积
     */
    public WktIntersectionResult intersection(String wgs84WKT1, String wgs84WKT2) {
        WktIntersectionResult result = new WktIntersectionResult();
        result.setWkt(config.EMPTYGEOM.toText());
        result.setMu(0.0);

        try {
            log.debug("开始计算两个WGS84几何图形的相交轮廓");

            // 1. 解析WKT字符串为几何对象
            Geometry geometry1 = toWgs84Geometry(wgs84WKT1);
            Geometry geometry2 = toWgs84Geometry(wgs84WKT2);

            log.debug("解析成功：几何1类型={}, 几何2类型={}", geometry1.getGeometryType(), geometry2.getGeometryType());

            // 2. 将WGS84几何图形转换为高斯投影坐标系以获得更精确的相交计算
            Geometry gaussGeometry1 = toGaussGeometry(geometry1);
            Geometry gaussGeometry2 = toGaussGeometry(geometry2);

            if (gaussGeometry1.isEmpty() || gaussGeometry2.isEmpty()) {
                log.warn("WGS84几何图形转换为高斯投影失败");
                return result;
            }

            log.debug("高斯投影转换成功：几何1类型={}, 几何2类型={}",
                    gaussGeometry1.getGeometryType(), gaussGeometry2.getGeometryType());

            // 3. 在高斯投影坐标系下执行相交操作（更精确）
            Geometry gaussIntersection = gaussGeometry1.intersection(gaussGeometry2);

            // 4. 检查是否有实际相交区域
            if (gaussIntersection == null || gaussIntersection.isEmpty()) {
                log.debug("两个几何图形没有相交区域");
                return result;
            }

            log.debug("相交成功，高斯投影下相交几何类型：{}，面积：{}平方米",
                    gaussIntersection.getGeometryType(), gaussIntersection.getArea());

            // 5. 将相交结果转换回WGS84坐标系
            Geometry wgs84Intersection = toWgs84Geometry(gaussIntersection);
            if (wgs84Intersection.isEmpty()) {
                log.warn("高斯投影相交结果转换回WGS84失败");
                return result;
            }

            // 6. 设置相交结果的WKT字符串
            result.setWkt(wgs84Intersection.toText());

            // 7. 在高斯投影坐标系下计算面积（更精确），然后转换为亩
            result.setMu(calcMu(wgs84Intersection));

            log.debug("相交轮廓计算完成：亩数={}亩（基于WGS84精确计算）", result.getMu());
        } catch (Exception e) {
            log.warn("相交计算失败：{}", e.getMessage());
        }

        return result;
    }

    /**
     * 计算WGS84坐标系下几何图形的球面面积（平方米）
     * <p>
     * 该方法根据几何图形类型自动选择合适的面积计算策略，支持单多边形和多多边形：
     * <ul>
     *   <li>对于单多边形(Polygon)，直接调用calculatePolygonSphericalArea计算</li>
     *   <li>对于多多边形(MultiPolygon)，累加所有子多边形的面积</li>
     * </ul>
     * 算法基于球面几何原理，精确考虑地球曲率，与Turf.js算法保持一致，确保跨平台结果一致性。
     * 适用于全球范围内的精确面积测量，不受区域大小限制。
     * </p>
     * 
     * @param wgs84Geometry WGS84坐标系下的几何图形，支持Polygon和MultiPolygon类型
     * @return 球面面积（平方米），始终为正值；对于不支持的几何类型返回0.0
     */
    public double calculateSphericalArea(Geometry wgs84Geometry) {
        double totalArea = 0.0;

        // 步骤1：根据几何图形类型选择计算策略
        if (wgs84Geometry instanceof Polygon) {
            // 单多边形处理
            totalArea = calculatePolygonSphericalArea((Polygon) wgs84Geometry);
            log.trace("单多边形面积: {}平方米", totalArea);
        } else if (wgs84Geometry instanceof MultiPolygon) {
            // 多多边形处理：累加所有子多边形面积
            MultiPolygon multiPolygon = (MultiPolygon) wgs84Geometry;
            for (int i = 0; i < multiPolygon.getNumGeometries(); i++) {
                Polygon polygon = (Polygon) multiPolygon.getGeometryN(i);
                double polyArea = calculatePolygonSphericalArea(polygon);
                totalArea += polyArea;
                log.trace("多边形{}面积: {}平方米", i, polyArea);
            }
            log.trace("MULTIPOLYGON总面积: {}平方米", totalArea);
        }

        // 步骤2：确保返回正值面积
        return Math.abs(totalArea);
    }

    /**
     * 计算WGS84坐标系下几何图形的面积（亩）
     * <p>
     * 该方法将WGS84坐标系下几何图形的球面面积转换为亩数，适用于土地测量、农业生产统计和资源评估。
     * 内部使用高精度球面面积计算算法，与Turf.js保持一致，结果进行四舍五入精确到4位小数。
     * 转换关系：1亩 = 2000/3 平方米 ≈ 666.6667平方米。
     * </p>
     * 
     * @param wgs84Geometry WGS84坐标系下的几何图形，支持Polygon和MultiPolygon类型
     * @return 几何图形的面积（亩），四舍五入保留4位小数；计算失败返回0.0
     */
    public double calcMu(Geometry wgs84Geometry) {
        try {
            // 步骤1：使用球面面积算法计算平方米面积
            double areaSqm = calculateSphericalArea(wgs84Geometry);

            // 步骤2：平方米转换为亩（1亩 = 2000/3平方米）
            // 四舍五入保留4位小数
            double mu = Math.round((areaSqm / (2000.0 / 3.0)) * 10000.0) / 10000.0;

            log.trace("计算几何图形面积（亩）: {}亩", mu);
            return mu;
        } catch (Exception e) {
            log.warn("WGS84几何图形计算亩数失败: {}", e.getMessage());
            return 0.0;
        }
    }

    /**
     * 计算WGS84坐标系下WKT字符串表示的几何图形面积（亩）
     * <p>
     * 该方法是calcMu(Geometry)的便捷重载版本，直接接受WKT字符串作为输入，
     * 内部自动解析WKT为几何对象并计算面积。适用于从数据库或文件中读取的WKT格式几何数据。
     * </p>
     * 
     * @param wgs84Wkt WGS84坐标系下的WKT字符串，如"POLYGON((...))"或"MULTIPOLYGON(((...)))"
     * @return 几何图形的面积（亩），四舍五入保留4位小数；解析或计算失败返回0.0
     * @see #calcMu(Geometry)
     */
    public double calcMu(String wgs84Wkt) {
        return calcMu(toWgs84Geometry(wgs84Wkt));
    }

    /**
     * 道路拆分处理 - 提取和优化田间作业轨迹
     * <p>
     * 该方法通过一系列处理步骤，从原始GPS轨迹中提取有效作业区域的几何轮廓：
     * <ol>
     *   <li>异常点位过滤（空时间、无效坐标等）</li>
     *   <li>冗余点优化（相同位置点去重）</li>
     *   <li>时间排序与坐标转换</li>
     *   <li>轨迹特征计算（距离、速度、方向角）</li>
     *   <li>高速移动识别与轨迹分段</li>
     *   <li>点抽稀与几何缓冲区生成</li>
     *   <li>多区域合并与边界优化</li>
     * </ol>
     * 主要用于农业机械作业轨迹分析、作业面积计算和作业质量评估。
     * </p>
     * 
     * @param wgs84Points WGS84坐标系下的轨迹点列表，每个点包含经纬度和时间信息
     * @param totalWidthM 作业幅宽（米），用于生成轨迹两侧的缓冲区宽度
     * @return 包含拆分结果的对象，包括作业轮廓几何图形、WKT格式、总面积和分段信息
     * @throws IllegalArgumentException 当参数验证失败时抛出
     * @see SplitRoadResult
     * @see OutlinePart
     */
    public SplitRoadResult splitRoad(List<TrackPoint> wgs84Points, double totalWidthM) {
        // 参数验证
        if (CollUtil.isEmpty(wgs84Points)) {
            throw new IllegalArgumentException("轨迹点列表不能为空");
        }
        if (totalWidthM < config.MIN_WORKING_WIDTH_M) {
            throw new IllegalArgumentException("幅宽（米）不能小于" + config.MIN_WORKING_WIDTH_M);
        }

        log.debug("参数：wgs84点数量={}, 总幅宽(米)={}", wgs84Points.size(), totalWidthM);

        // 计算半宽（用于生成缓冲区）
        double halfWidthM = totalWidthM / 2.0;

        // 初始化结果对象
        SplitRoadResult result = new SplitRoadResult();
        result.setTotalWidthM(totalWidthM);
        result.setOutline(config.EMPTYGEOM);
        result.setWkt(config.EMPTYGEOM.toText());

        // 步骤1：过滤异常点位信息
        wgs84Points = wgs84Points.stream()
                .filter(p -> {
                    // 1. 时间不能为空
                    if (p.getTime() == null) {
                        log.warn("轨迹点时间为空，此点抛弃");
                        return false;
                    }
                    // 2. 经纬度不能为0（无效坐标）
                    if (p.getLon() == 0.0 && p.getLat() == 0.0) {
                        log.warn("定位时间: {} 轨迹点经纬度为0，此点抛弃", p.getTime());
                        return false;
                    }
                    // 3. 经纬度必须在合理范围内
                    if (p.getLon() < -180.0 || p.getLon() > 180.0 || p.getLat() < -90.0 || p.getLat() > 90.0) {
                        log.warn("定位时间: {} 轨迹点经纬度超出范围：[{},{}] 此点抛弃", p.getTime(), p.getLon(), p.getLat());
                        return false;
                    }
                    return true;
                })
                .collect(Collectors.toList());
        log.debug("过滤异常点位后轨迹点数量={}", wgs84Points.size());

        // 步骤2：过滤相邻经纬度完全一致的点位（去重，减少计算量）
        log.debug("过滤相邻经纬度完全一致的点位");
        if (wgs84Points.size() > 1) {
            List<TrackPoint> filteredPoints = new ArrayList<>();
            filteredPoints.add(wgs84Points.get(0)); // 保留第一个点

            for (int i = 1; i < wgs84Points.size(); i++) {
                TrackPoint current = wgs84Points.get(i);
                TrackPoint previous = filteredPoints.get(filteredPoints.size() - 1);

                // 比较经纬度是否完全一致（使用极小值进行比较，避免浮点精度问题）
                if (Math.abs(current.getLon() - previous.getLon()) < 1e-10 &&
                        Math.abs(current.getLat() - previous.getLat()) < 1e-10) {
                    // 记录被过滤掉的点的定位时间
                    log.trace("过滤掉经纬度与前一点相同的点位，定位时间: {}", current.getTime());
                } else {
                    filteredPoints.add(current);
                }
            }

            wgs84Points = filteredPoints;
        }
        log.debug("过滤后轨迹点数量={}", wgs84Points.size());

        // 步骤3：按定位时间升序排序，确保轨迹时序正确
        wgs84Points.sort(Comparator.comparing(TrackPoint::getTime));

        // 获取上报时间间隔分布
        Map<Integer, Integer> intervalDistribution = new HashMap<>();
        log.trace("统计轨迹点时间间隔分布");
        for (int i = 1; i < wgs84Points.size(); i++) {
            TrackPoint prevPoint = wgs84Points.get(i - 1);
            TrackPoint currPoint = wgs84Points.get(i);

            // 计算时间间隔（秒）- LocalDateTime使用Duration
            Duration duration = Duration.between(prevPoint.getTime(), currPoint.getTime());
            int timeDiffSeconds = (int) duration.getSeconds();
            if (timeDiffSeconds < 12) {
                // 统计每个间隔的出现次数
                intervalDistribution.put(timeDiffSeconds, intervalDistribution.getOrDefault(timeDiffSeconds, 0) + 1);
            }
        }
        // 获取最小有效时间间隔
        int minEffectiveInterval = 1; // 默认值
        if (!intervalDistribution.isEmpty()) {
            // 找到点数最多的时间间隔，如果点数相同则选择时间间隔更小的
            minEffectiveInterval = intervalDistribution.entrySet().stream()
                    .max((e1, e2) -> {
                        int countCompare = Integer.compare(e1.getValue(), e2.getValue());
                        if (countCompare != 0) {
                            return countCompare; // 点数多的优先
                        }
                        return Integer.compare(e2.getKey(), e1.getKey()); // 点数相同时，时间间隔小的优先（降序比较）
                    })
                    .map(Map.Entry::getKey)
                    .orElse(1);
        }
        log.info("最小有效时间间隔：{}秒", minEffectiveInterval);

        // 创建final变量供lambda表达式使用
        final int finalMinEffectiveInterval = minEffectiveInterval;

        // 步骤4：转换到高斯投影平面坐标，便于后续距离和几何计算
        List<TrackPoint> gaussPoints = toGaussPointList(wgs84Points);

        // 步骤5：填充轨迹点的其他信息：距离、速度、方向角
        for (int i = 1; i < gaussPoints.size(); i++) {
            TrackPoint prevPoint = gaussPoints.get(i - 1);
            TrackPoint currPoint = gaussPoints.get(i);

            // 计算两点间的距离（米）- 高斯投影下直接使用欧几里得距离
            double deltaX = currPoint.getLon() - prevPoint.getLon();
            double deltaY = currPoint.getLat() - prevPoint.getLat();
            double distance = Math.sqrt(deltaX * deltaX + deltaY * deltaY);
            currPoint.setDistance(distance);

            // 计算两点间的时间差（秒）
            long timeDiffSeconds = 0;
            if (prevPoint.getTime() != null && currPoint.getTime() != null) {
                timeDiffSeconds = Duration.between(prevPoint.getTime(), currPoint.getTime()).getSeconds();
            }

            // 计算速度（米/秒）
            double speed = 0;
            if (timeDiffSeconds > 0) {
                speed = distance / timeDiffSeconds;
            }
            currPoint.setSpeed(speed);

            // 计算方向角（度）
            double direction = 0;
            if (distance > 0) {
                // 使用atan2计算弧度方向，然后转换为角度
                direction = Math.toDegrees(Math.atan2(deltaY, deltaX));
                // 确保方向角在0-360度范围内
                if (direction < 0) {
                    direction += 360;
                }
            }
            currPoint.setDirection(direction);
        }
        // 循环打印所有轨迹点的信息（调试代码，已注释）
        /* for (TrackPoint point : gaussPoints) {
            log.debug("定位时间：{} 轨迹点：[{},{}] 距离：{} 速度：{} 方向：{}",
                    point.getTime(), point.getLon(), point.getLat(), point.getDistance(), point.getSpeed(),
                    point.getDirection());
        } */

        // 步骤6：根据策略进行轨迹分段（识别作业与非作业状态）
        log.trace("切分成多段轨迹");
        List<List<TrackPoint>> gaussTrackSegments = new ArrayList<>();
        List<TrackPoint> currentSegment = new ArrayList<>();
        currentSegment.add(gaussPoints.get(0)); // 添加第一个点到当前段

        for (int i = 1; i < gaussPoints.size(); i++) {
            TrackPoint point = gaussPoints.get(i);
            TrackPoint prevPoint = gaussPoints.get(i - 1);

            // 计算时间间隔（秒）
            long timeInterval = Duration.between(prevPoint.getTime(), point.getTime()).getSeconds();

            // 高速通常表示转移而非作业，或者时间间隔过大也表示可能有中断
            // km/h 转换为 m/s 的公式：kmh / 3.6 = m/s
            // 增加条件：如果时间间隔超过minEffectiveInterval的3倍，也切割新段
            double max_work_speed;
            if (minEffectiveInterval < 5) {
                max_work_speed = config.MAX_WORK_SPEED_1s / 3.6;
            } else {
                max_work_speed = config.MAX_WORK_SPEED_10s / 3.6;
            }
            if (point.getSpeed() > max_work_speed
                    || timeInterval > minEffectiveInterval * 3) {
                // 如果当前段不为空且有多个点，则将当前段添加到结果中
                if (currentSegment.size() > 1) {
                    gaussTrackSegments.add(new ArrayList<>(currentSegment));
                    log.trace("创建新轨迹段，上一段包含{}个点", currentSegment.size());
                }
                // 开始新的段
                currentSegment.clear();
            }

            // 将当前点添加到当前段
            currentSegment.add(point);
        }

        // 添加最后一个轨迹段（如果不为空且有多个点）
        if (currentSegment.size() > 1) {
            gaussTrackSegments.add(currentSegment);
            log.trace("添加最后一个轨迹段，包含{}个点", currentSegment.size());
        }
        log.debug("轨迹切割完成，共得到 {} 段轨迹", gaussTrackSegments.size());

        // 步骤7：处理每个轨迹段，生成几何轮廓
        List<Geometry> gaussGeometries = new ArrayList<>();
        for (List<TrackPoint> segment : gaussTrackSegments) {
            log.trace("处理轨迹段，原始点数: {}", segment.size());

            long startTime = System.currentTimeMillis();

            // 7.1 点抽稀优化 - 使用Douglas-Peucker算法减少点数但保持形状特征
            List<TrackPoint> simplifiedPoints = simplifyTrackPoints(segment, 0.1); // 0.1米容差
            log.trace("点抽稀后点数: {}, 耗时: {}ms", simplifiedPoints.size(), System.currentTimeMillis() - startTime);

            // 7.2 根据点数选择处理策略 - 大规模轨迹采用分块处理
            Geometry gaussGeometry;
            if (simplifiedPoints.size() > 1000) {
                // 当点数超过1000时，进行分块处理以避免内存溢出和提高性能
                gaussGeometry = processLargeSegmentInChunks(simplifiedPoints, halfWidthM);
            } else {
                // 点数较少时直接处理：创建线串，然后生成缓冲区（表示作业范围）
                Coordinate[] coordinates = simplifiedPoints.stream()
                        .map(point -> new Coordinate(point.getLon(), point.getLat()))
                        .toArray(Coordinate[]::new);
                LineString lineString = config.geometryFactory.createLineString(coordinates);
                // 创建宽度为halfWidthM的缓冲区，表示作业区域
                gaussGeometry = lineString.buffer(halfWidthM);
            }
            gaussGeometries.add(gaussGeometry);
        }

        // 将gaussGeometries列表转换为Geometry[]数组
        Geometry unionGaussGeometry = config.geometryFactory
                .createGeometryCollection(gaussGeometries.toArray(new Geometry[0]))
                .union()
                .buffer(0.03).buffer(-0.03);

        List<OutlinePart> outlineParts = new ArrayList<>();
        if (unionGaussGeometry instanceof Polygon) {
            Geometry wgs84Geometry = toWgs84Geometry(unionGaussGeometry);
            // 7.4 性能优化：获取几何图形的边界框，用于快速空间过滤
            Envelope geometryEnvelope = wgs84Geometry.getEnvelopeInternal();
            // 使用PreparedGeometry预优化，大幅提升空间判断性能（10-100倍）
            PreparedGeometry preparedWgs84Geometry = PreparedGeometryFactory.prepare(wgs84Geometry);
            List<TrackPoint> wgs84PointsSegment = wgs84Points.stream()
                    .filter(point -> {
                        try {
                            // 第一步：边界框快速过滤 - 使用简单的数值比较，性能极高
                            if (!geometryEnvelope.contains(point.getLon(), point.getLat())) {
                                return false;
                            }
                            // 第二步：使用PreparedGeometry进行精确空间判断
                            return preparedWgs84Geometry.contains(config.geometryFactory.createPoint(
                                    new Coordinate(point.getLon(), point.getLat())));
                        } catch (Exception e) {
                            log.trace("点位空间判断失败：经度{} 纬度{} 错误：{}",
                                    point.getLon(), point.getLat(), e.getMessage());
                            return false;
                        }
                    })
                    .collect(Collectors.toList());
            log.trace("点位空间判断完成，原始点位数：{}，筛选后点位数：{}", wgs84Points.size(), wgs84PointsSegment.size());

            OutlinePart outlinePart = new OutlinePart();
            outlinePart.setTotalWidthM(totalWidthM);
            outlinePart.setOutline(unionGaussGeometry);
            outlinePart.setWkt(wgs84Geometry.toText());
            outlinePart.setMu(calcMu(wgs84Geometry));
            outlinePart.setTrackPoints(wgs84PointsSegment);
            outlinePart.setStartTime(wgs84PointsSegment.get(0).getTime());
            outlinePart.setEndTime(wgs84PointsSegment.get(wgs84PointsSegment.size() - 1).getTime());
            outlineParts.add(outlinePart);
        } else if (unionGaussGeometry instanceof MultiPolygon) {
            for (int i = 0; i < unionGaussGeometry.getNumGeometries(); i++) {
                Geometry partGaussGeometry = unionGaussGeometry.getGeometryN(i);
                Geometry wgs84GeometryPart = toWgs84Geometry(partGaussGeometry);
                Envelope geometryEnvelopePart = wgs84GeometryPart.getEnvelopeInternal();
                PreparedGeometry preparedWgs84GeometryPart = PreparedGeometryFactory.prepare(wgs84GeometryPart);
                List<TrackPoint> wgs84PointsSegmentPart = wgs84Points.stream()
                        .filter(point -> {
                            try {
                                // 第一步：边界框快速过滤 - 使用简单的数值比较，性能极高
                                if (!geometryEnvelopePart.contains(point.getLon(), point.getLat())) {
                                    return false;
                                }
                                // 第二步：使用PreparedGeometry进行精确空间判断
                                return preparedWgs84GeometryPart.contains(config.geometryFactory.createPoint(
                                        new Coordinate(point.getLon(), point.getLat())));
                            } catch (Exception e) {
                                log.trace("点位空间判断失败：经度{} 纬度{} 错误：{}",
                                        point.getLon(), point.getLat(), e.getMessage());
                                return false;
                            }
                        })
                        .collect(Collectors.toList());
                log.trace("点位空间判断完成，原始点位数：{}，筛选后点位数：{}", wgs84Points.size(), wgs84PointsSegmentPart.size());

                OutlinePart outlinePart = new OutlinePart();
                outlinePart.setTotalWidthM(totalWidthM);
                outlinePart.setOutline(partGaussGeometry);
                outlinePart.setWkt(wgs84GeometryPart.toText());
                outlinePart.setMu(calcMu(wgs84GeometryPart));
                outlinePart.setTrackPoints(wgs84PointsSegmentPart);
                outlinePart.setStartTime(wgs84PointsSegmentPart.get(0).getTime());
                outlinePart.setEndTime(wgs84PointsSegmentPart.get(wgs84PointsSegmentPart.size() - 1).getTime());
                outlineParts.add(outlinePart);
            }
        }

        // 按outlineParts里的亩数降序排序，然后只保留亩数最高的10个多边形
        log.trace("按亩数降序排序");
        outlineParts.sort(Comparator.comparingDouble(OutlinePart::getMu).reversed());

        if (outlineParts.size() > config.MAX_GEOMETRY) {
            log.debug("截取前有 {} 个多边形，只保留亩数最大的{}个", outlineParts.size(), config.MAX_GEOMETRY);
            outlineParts = outlineParts.subList(0, Math.min(config.MAX_GEOMETRY, outlineParts.size()));
        }

        if (outlineParts.size() > 1) {
            if (minEffectiveInterval < 5) {
                outlineParts = outlineParts.stream()
                        .filter(part -> part.getTrackPoints().size() >= config.MIN_WORK_POINTS_1s)
                        .collect(Collectors.toList());
            } else {
                outlineParts = outlineParts.stream()
                        .filter(part -> part.getTrackPoints().size() >= config.MIN_WORK_POINTS_10s)
                        .collect(Collectors.toList());
            }
            log.debug("过滤最小点位，剩余 {} 个多边形", outlineParts.size());
        }

        // 过滤掉亩数低于MIN_MU的多边形
        if (outlineParts.size() > 1) {
            OutlinePart bak = outlineParts.get(0);
            outlineParts = outlineParts.stream()
                    .filter(part -> part
                            .getMu() >= (finalMinEffectiveInterval < 5 ? config.MIN_MU_1s : config.MIN_MU_10s))
                    .collect(Collectors.toList());
            if (outlineParts.isEmpty()) {
                outlineParts.add(bak);
            }
            log.debug("过滤最小亩数，剩余 {} 个多边形", outlineParts.size());
        }

        // 合并outlineParts中的所有outline几何图形
        unionGaussGeometry = config.geometryFactory
                .createGeometryCollection(outlineParts.stream()
                        .map(OutlinePart::getOutline)
                        .toArray(Geometry[]::new))
                .union();

        // 步骤9：计算总面积（亩）
        double totalMu = outlineParts != null ? outlineParts.stream()
                .filter(Objects::nonNull)
                .mapToDouble(OutlinePart::getMu)
                .sum() : 0.0;

        // 步骤10：设置最终结果
        result.setOutline(unionGaussGeometry);
        result.setWkt(toWgs84Geometry(unionGaussGeometry).toText());
        result.setParts(outlineParts);
        result.setMu(Math.round(totalMu * 10000.0) / 10000.0); // 四舍五入保留4位小数

        log.debug("最终生成区块数量：{}块，总亩数：{}亩，最终几何类型：{}",
                unionGaussGeometry.getNumGeometries(), result.getMu(), unionGaussGeometry.getGeometryType());

        return result;
    }

}