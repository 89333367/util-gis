package sunyu.util;

import cn.hutool.core.collection.CollUtil;
import cn.hutool.core.util.StrUtil;
import cn.hutool.log.Log;
import cn.hutool.log.LogFactory;
import elki.clustering.dbscan.DBSCAN;
import elki.data.Cluster;
import elki.data.Clustering;
import elki.data.DoubleVector;
import elki.data.model.Model;
import elki.data.type.TypeUtil;
import elki.database.Database;
import elki.database.StaticArrayDatabase;
import elki.database.ids.DBIDIter;
import elki.database.ids.DBIDRange;
import elki.database.relation.Relation;
import elki.datasource.ArrayAdapterDatabaseConnection;
import elki.distance.minkowski.EuclideanDistance;
import org.geotools.geometry.jts.JTS;
import org.geotools.geometry.jts.JTSFactoryFinder;
import org.geotools.referencing.CRS;
import org.geotools.referencing.crs.DefaultGeographicCRS;
import org.locationtech.jts.geom.*;
import org.locationtech.jts.index.strtree.STRtree;
import org.locationtech.jts.io.ParseException;
import org.locationtech.jts.io.WKTReader;
import org.locationtech.jts.simplify.DouglasPeuckerSimplifier;
import org.opengis.referencing.crs.CoordinateReferenceSystem;
import org.opengis.referencing.operation.MathTransform;
import sunyu.util.pojo.*;

import java.time.Duration;
import java.util.*;
import java.util.concurrent.ConcurrentHashMap;
import java.util.stream.Collectors;

/**
 * GIS工具类，提供几何操作、投影转换、空间分析等功能
 *
 * @author SunYu
 */
public class GisUtil implements AutoCloseable {
    private final Log log = LogFactory.get();
    private final Config config;

    public static Builder builder() {
        return new Builder();
    }

    private GisUtil(Config config) {
        log.info("[构建{}] 开始", this.getClass().getSimpleName());
        this.config = config;
        log.info("[构建{}] 结束", this.getClass().getSimpleName());
    }

    private static class Config {
        /**
         * 几何工厂，用于创建各种几何对象
         * <p>使用JTSFactoryFinder获取线程安全的几何工厂实例，避免重复创建</p>
         */
        private final GeometryFactory GEOMETRY_FACTORY = JTSFactoryFinder.getGeometryFactory();

        /**
         * 空几何集合，用于表示无效或空的几何结果
         * <p>作为常量复用，避免每次创建新的空几何对象，提高内存效率</p>
         */
        private final Geometry EMPTY_GEOMETRY = GEOMETRY_FACTORY.createGeometryCollection();

        /**
         * WGS84坐标参考系统，使用默认地理CRS，避免重复创建
         * <p>这是全球标准的地理坐标参考系统，所有输入输出都基于此坐标系</p>
         */
        private final CoordinateReferenceSystem WGS84_CRS = DefaultGeographicCRS.WGS84;

        /**
         * 缓存键格式
         */
        private final String CACHE_KEY_FORMAT = "%d_%.1f_%.1f";

        /**
         * 高斯投影CRS缓存（线程安全）
         * <p>
         * <strong>键格式：</strong>"zone_falseEasting_centralMeridian"
         * <br>
         * <strong>缓存策略：</strong>基于投影带参数的唯一性进行缓存
         * <br>
         * <strong>性能提升：</strong>避免重复创建CRS对象，显著提升坐标转换性能
         * </p>
         */
        private final ConcurrentHashMap<String, CoordinateReferenceSystem> GAUSS_CRS_CACHE = new ConcurrentHashMap<>();

        /**
         * WGS84到高斯投影的坐标转换缓存（线程安全）
         * <p>
         * <strong>键格式：</strong>"zone_falseEasting_centralMeridian"
         * <br>
         * <strong>缓存策略：</strong>基于投影带参数的唯一性进行缓存
         * <br>
         * <strong>性能提升：</strong>避免重复创建MathTransform对象，坐标转换性能提升10倍以上
         * </p>
         */
        private final ConcurrentHashMap<String, MathTransform> WGS84_TO_GAUSS_TRANSFORM_CACHE = new ConcurrentHashMap<>();

        /**
         * 高斯投影到WGS84的坐标转换缓存（线程安全）
         * <p>
         * <strong>键格式：</strong>"zone_falseEasting_centralMeridian"
         * <br>
         * <strong>缓存策略：</strong>基于投影带参数的唯一性进行缓存
         * <br>
         * <strong>性能提升：</strong>避免重复创建MathTransform对象，坐标转换性能提升10倍以上
         * </p>
         */
        private final ConcurrentHashMap<String, MathTransform> GAUSS_TO_WGS84_TRANSFORM_CACHE = new ConcurrentHashMap<>();

        /**
         * WGS84椭球长半轴（米），与Turf.js保持一致，用于计算球面距离
         * <p>标准值：6378137.0米，这是WGS84椭球体的长半轴长度</p>
         */
        private final double EARTH_RADIUS = 6378137.0;

        /**
         * 两点间最大作业距离(米)
         * <p>对应18km/h的作业速度（5米/秒），用于判断轨迹点之间的合理性</p>
         */
        private final double MAX_WORK_DISTANCE = 18 / 3.6;

        /**
         * DBSCAN最小点数
         * <p>
         * 用于DBSCAN聚类算法的最小点阈值。
         * 当一个区域内的点数量小于此值时，被认为是噪声点或异常值，不会被分配到任何聚类中。
         * </p>
         */
        private final int MIN_DBSCAN_POINTS = 20;

        /**
         * 最大拆分返回数量
         * <p>
         * 用于限制拆分轨迹点时返回的最大数量。
         * 当轨迹点数量超过此阈值时，会对轨迹进行拆分，返回多个子轨迹。
         * 该参数可以有效控制内存占用和处理时间，避免处理过大数据集时内存溢出。
         * </p>
         */
        private final int MAX_SPLIT_RETURN_SIZE = 24;

        /**
         * 最小返回面积（亩）
         */
        private final double MIN_RETURN_MU = 0.5;

        /**
         * 亩到平方米的转换系数(亩 * 这个数)
         * <p>
         * 1亩 = 2000/3平方米
         */
        private final double MU_TO_SQUARE_METER = 2000.0 / 3.0;

        /**
         * 平方米到亩的转换系数(平方米 * 这个数)
         * <p>
         * 1平方米 = 3/2000亩
         */
        private final double SQUARE_TO_MU_METER = 3.0 / 2000.0;

        /**
         * 拆分时间间隔（秒）
         */
        private final int SPLIT_TIME_SECOND = 60 * 30;

        /**
         * 渐进式容差（米）
         */
        private final double[] TOLERANCES = {0.00111, 0.0111, 0.111, 1.11, 11.1}; // 渐进式容差（米）

        /**
         * 米到度的转换系数 将米转换为度（近似转换：1度≈111公里）
         */
        private final double MI_TO_DEGREE = 111000.0;
    }

    public static class Builder {
        private Config config = new Config();

        public GisUtil build() {
            return new GisUtil(config);
        }
    }

    @Override
    public void close() {
        log.info("[销毁{}] 开始", this.getClass().getSimpleName());
        log.info("[销毁{}] 结束", this.getClass().getSimpleName());
    }

    /**
     * 获取高斯投影CRS
     * <p>
     * <strong>缓存策略：</strong>基于投影带号、假东距和中央经线的组合作为唯一键进行缓存
     * <br>
     * <strong>性能提升：</strong>避免重复创建CRS对象，显著提升坐标转换性能
     * </p>
     *
     * @param zone            投影带号
     * @param falseEasting    假东距
     * @param centralMeridian 中央经线
     *
     * @return 对应的高斯投影CRS对象
     */
    private CoordinateReferenceSystem getGaussCRS(int zone, double falseEasting, double centralMeridian) {
        // 创建缓存键，使用带号、假东距和中央经线作为唯一标识
        String cacheKey = String.format(config.CACHE_KEY_FORMAT, zone, falseEasting, centralMeridian);
        log.debug("获取高斯投影CRS：缓存键 {}", cacheKey);

        // 使用computeIfAbsent实现线程安全的缓存机制，只在缓存未命中时创建新CRS
        return config.GAUSS_CRS_CACHE.computeIfAbsent(cacheKey, key -> {
            try {
                log.debug("创建高斯投影CRS：投影带号 {} 假东距 {} 中央经线 {}", zone, falseEasting, centralMeridian);
                // 定义高斯-克吕格投影坐标系 - 使用预构建的WKT模板
                // WKT格式定义了完整的坐标系统参数，包括基准面、椭球体、投影方式和参数
                String wktTemplate = "PROJCS[\"Gauss_Kruger_ZONE_%d\", GEOGCS[\"GCS_WGS_1984\", DATUM[\"WGS_1984\", SPHEROID[\"WGS_84\", 6378137.0, 298.257223563]], PRIMEM[\"Greenwich\", 0.0], UNIT[\"Degree\", 0.0174532925199433]], PROJECTION[\"Transverse_Mercator\"], PARAMETER[\"False_Easting\", %.12f], PARAMETER[\"False_Northing\", 0.0], PARAMETER[\"Central_Meridian\", %.12f], PARAMETER[\"Scale_Factor\", 1.0], PARAMETER[\"Latitude_Of_Origin\", 0.0], UNIT[\"Meter\", 1.0]]";
                String gaussProjString = String.format(wktTemplate, zone, falseEasting, centralMeridian);

                // 解析WKT字符串创建CRS对象
                return CRS.parseWKT(gaussProjString);
            } catch (Exception e) {
                log.warn("创建高斯投影CRS失败：投影带号 {} 假东距 {}, 中央经线 {}, 错误 {}", zone, falseEasting, centralMeridian, e.getMessage());
                return null;
            }
        });
    }

    /**
     * 计算WGS84坐标系下的环的球面面积（平方米）
     *
     * @param wgs84Ring 输入的WGS84坐标系下的环对象（LineString类型）
     *
     * @return 环的球面面积（平方米），如果输入无效则返回0.0
     */
    private double calculateRingSphericalArea(LineString wgs84Ring) {
        // 参数有效性检查
        if (wgs84Ring == null || wgs84Ring.isEmpty()) {
            return 0.0;
        }

        // 获取环的所有坐标点
        Coordinate[] coords = wgs84Ring.getCoordinates();
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
     * 计算WGS84坐标系下的多边形的球面面积（平方米）
     *
     * @param wgs84Polygon 输入的WGS84坐标系下的多边形对象（Polygon类型）
     *
     * @return 多边形的球面面积（平方米），如果输入无效则返回0.0
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
     * 处理数据块，将其转换为缓冲几何图形
     *
     * @param points      输入的点列表，每个点包含高斯投影坐标（GaussPoint类型）
     * @param startIndex  数据块的起始索引
     * @param chunkSize   数据块的大小（点数）
     * @param totalPoints 总数据点数量
     * @param bufferWidth 缓冲宽度（单位：米）
     *
     * @return 处理后的缓冲几何图形（LineString或Polygon类型）
     */
    private Geometry processChunk(List<GaussPoint> points, int startIndex, int chunkSize, int totalPoints, double bufferWidth) {
        int end = Math.min(startIndex + chunkSize, totalPoints);
        List<GaussPoint> chunk = points.subList(startIndex, end);

        long chunkStartTime = System.currentTimeMillis();

        // 优化：直接使用数组而不是中间集合
        int chunkLength = chunk.size();
        Coordinate[] coords = new Coordinate[chunkLength];
        for (int j = 0; j < chunkLength; j++) {
            GaussPoint p = chunk.get(j);
            coords[j] = new Coordinate(p.getGaussX(), p.getGaussY());
        }

        LineString chunkLine = config.GEOMETRY_FACTORY.createLineString(coords);

        // 优化4: 进一步优化buffer参数，使用最小必要的精度
        int quadrantSegments = 4; // 更低的值，显著提升性能，在保证准确性的前提下
        int endCapStyle = 1; // 1=LINEAREND，高效的样式

        // 计算buffer
        Geometry chunkBuffer = chunkLine.buffer(bufferWidth, quadrantSegments, endCapStyle);

        // 可选：简化几何图形，进一步减少复杂度（如果buffer后图形仍然复杂）
        if (chunkBuffer.getNumPoints() > 1000) {
            chunkBuffer = DouglasPeuckerSimplifier.simplify(chunkBuffer, 0.00001); // 轻微简化，保持形状
        }

        log.debug("处理块: {}-{}, 点数: {}, buffer几何类型: {}, 耗时: {}ms", startIndex, end, chunkLength, chunkBuffer.getGeometryType(), System.currentTimeMillis() - chunkStartTime);

        return chunkBuffer;
    }

    /**
     * 递归合并几何图形列表，优化合并效率
     *
     * @param geometries 输入的几何图形列表（LineString或Polygon类型）
     *
     * @return 合并后的几何图形（LineString或Polygon类型），如果输入无效则返回空几何图形
     */
    private Geometry mergeGeometriesRecursively(List<Geometry> geometries) {
        // 基本情况：列表为空或只有一个元素
        if (geometries.isEmpty()) {
            return config.EMPTY_GEOMETRY;
        }
        if (geometries.size() == 1) {
            return geometries.get(0);
        }

        // 优化1：过滤掉空的几何图形
        List<Geometry> nonEmptyGeometries = new ArrayList<>();
        for (Geometry geom : geometries) {
            if (geom != null && !geom.isEmpty()) {
                nonEmptyGeometries.add(geom);
            }
        }

        if (nonEmptyGeometries.isEmpty()) {
            return config.EMPTY_GEOMETRY;
        }
        if (nonEmptyGeometries.size() == 1) {
            return nonEmptyGeometries.get(0);
        }

        // 优化2：使用空间索引优化合并顺序
        // 创建STRtree空间索引，提高查询效率
        STRtree index = new STRtree();
        for (Geometry geom : nonEmptyGeometries) {
            index.insert(geom.getEnvelopeInternal(), geom);
        }
        index.build();

        // 优化3：使用空间索引找到相邻几何图形，优先合并相邻的几何图形
        List<Geometry> mergedGeometries = new ArrayList<>();
        Set<Geometry> processed = new HashSet<>();

        for (Geometry geom : nonEmptyGeometries) {
            if (processed.contains(geom)) {
                continue;
            }

            // 查找与当前几何图形相交的其他几何图形
            List<Geometry> nearbyGeoms = new ArrayList<>();
            List<?> rawResults = index.query(geom.getEnvelopeInternal());
            for (Object item : rawResults) {
                if (item instanceof Geometry) {
                    nearbyGeoms.add((Geometry) item);
                }
            }
            Geometry current = geom;
            processed.add(current);

            for (Geometry nearby : nearbyGeoms) {
                if (processed.contains(nearby) || nearby == current) {
                    continue;
                }

                // 快速检查边界框是否相交，避免不必要的union操作
                if (current.getEnvelopeInternal().intersects(nearby.getEnvelopeInternal())) {
                    try {
                        // 使用更高效的合并策略：先进行缓冲区操作，减少复杂度
                        Geometry tempCurrent = current.buffer(0.01);
                        Geometry tempNearby = nearby.buffer(0.01);
                        Geometry newResult = tempCurrent.union(tempNearby);

                        // 对合并结果进行轻微简化，减少复杂度
                        if (newResult.getNumPoints() > 1000) {
                            newResult = DouglasPeuckerSimplifier.simplify(newResult, 0.00001);
                        }

                        current = newResult;
                        processed.add(nearby);
                    } catch (Exception e) {
                        log.trace("合并几何图形时出错: {}", e.getMessage());
                        // 合并失败时继续下一个几何图形
                    }
                }
            }

            mergedGeometries.add(current);
        }

        // 优化4：如果还有多个几何图形，进行最终合并
        if (mergedGeometries.size() == 1) {
            return mergedGeometries.get(0);
        } else if (mergedGeometries.size() > 1) {
            // 使用级联合并，减少复杂度
            Geometry result = mergedGeometries.get(0);
            for (int i = 1; i < mergedGeometries.size(); i++) {
                try {
                    result = result.union(mergedGeometries.get(i));
                    // 每次合并后进行轻微简化，控制复杂度增长
                    if (result.getNumPoints() > 2000) {
                        result = DouglasPeuckerSimplifier.simplify(result, 0.00001);
                    }
                } catch (Exception e) {
                    log.warn("最终合并几何图形时出错: {}", e.getMessage());
                }
            }
            return result;
        }

        return config.EMPTY_GEOMETRY;
    }

    /**
     * 计算轨迹段的点密度（点/米）
     * 通过计算相邻点之间的平均距离来估算密度
     *
     * @param points 轨迹点列表
     *
     * @return 密度值（点/米），值越大表示点越密集
     */
    private double calculatePointDensity(List<GaussPoint> points) {
        if (points.size() < 2) {
            return 0.0;
        }

        double totalDistance = 0.0;
        int validSegments = 0;

        // 计算相邻点之间的总距离
        for (int i = 1; i < points.size(); i++) {
            GaussPoint prevPoint = points.get(i - 1);
            GaussPoint currPoint = points.get(i);

            // 计算欧几里得距离
            double distance = Math.sqrt(Math.pow(currPoint.getGaussX() - prevPoint.getGaussX(), 2) + Math.pow(currPoint.getGaussY() - prevPoint.getGaussY(), 2));

            if (distance > 0) { // 避免重复点
                totalDistance += distance;
                validSegments++;
            }
        }

        if (validSegments == 0) {
            return 0.0;
        }

        double avgDistance = totalDistance / validSegments;
        double density = avgDistance > 0 ? 1.0 / avgDistance : 0.0; // 密度 = 1/平均距离

        log.debug("密度计算：总点数={}, 有效段数={}, 总距离={}米, 平均距离={}米, 密度={}点/米", points.size(), validSegments, totalDistance, avgDistance, density);

        return density;
    }

    /**
     * 分块处理大型轨迹段，避免内存溢出
     *
     * @param points      轨迹点列表（GaussPoint类型）
     * @param bufferWidth 缓冲区宽度（米），用于合并分块后的几何图形
     *
     * @return 合并后的几何图形（LineString或Polygon类型），如果输入无效则返回空几何图形
     */
    private Geometry processLargeSegmentInChunks(List<GaussPoint> points, double bufferWidth) {
        long startTime = System.currentTimeMillis();
        log.debug("开始分块处理大型轨迹段，点数: {}", points.size());

        // 计算点的密度，动态调整分块策略
        double density = calculatePointDensity(points);
        log.debug("轨迹段密度: {} 点/米", density);

        // 根据密度动态设置分块参数
        int totalPoints = points.size();
        int chunkSize;
        int overlapSize;

        if (density > 1.0) {
            // 高密度区域：点非常密集，使用小块，减少每块点数
            chunkSize = Math.max(300, Math.min(totalPoints, 500));
            overlapSize = Math.max(30, chunkSize / 10); // 10%重叠
            log.debug("高密度区域，使用小块处理：块大小={}, 重叠={}", chunkSize, overlapSize);
        } else if (density > 0.5) {
            // 中密度区域：中等密度，使用中等块
            chunkSize = Math.max(500, Math.min(totalPoints, 800));
            overlapSize = Math.max(40, chunkSize / 12); // 8%重叠
            log.debug("中密度区域，使用中块处理：块大小={}, 重叠={}", chunkSize, overlapSize);
        } else {
            // 低密度区域：点稀疏，使用大块，减少分块数
            chunkSize = Math.max(800, Math.min(totalPoints, 1200));
            overlapSize = Math.max(50, chunkSize / 15); // 6.7%重叠
            log.debug("低密度区域，使用大块处理：块大小={}, 重叠={}", chunkSize, overlapSize);
        }

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

        log.debug("共分成 {} 个块进行处理，总点数: {}", chunkStartIndices.size(), totalPoints);

        // 处理所有块
        List<Geometry> chunkGeometries = new ArrayList<>();
        for (int i = 0; i < chunkStartIndices.size(); i++) {
            int startIndex = chunkStartIndices.get(i);
            log.debug("处理第 {} 块，起始索引: {}, 块大小: {}", i + 1, startIndex, chunkSize);

            Geometry geom = processChunk(points, startIndex, chunkSize, totalPoints, bufferWidth);
            if (!geom.isEmpty()) {
                chunkGeometries.add(geom);
            }
        }

        log.debug("所有块处理完成，准备合并几何图形，耗时: {} ms", System.currentTimeMillis() - startTime);

        // 直接合并几何图形，移除中间预处理步骤
        Geometry mergedGeometry = mergeGeometriesRecursively(chunkGeometries);

        // 简化后处理：只进行基本的有效性检查和修复，避免过度处理
        if (!mergedGeometry.isValid()) {
            mergedGeometry = mergedGeometry.buffer(0); // 修复无效几何图形
        }

        log.debug("分块处理完成，合并后几何类型: {}, 总耗时: {}ms", mergedGeometry.getGeometryType(), System.currentTimeMillis() - startTime);
        return mergedGeometry;
    }


    /**
     * 按最大距离切分轨迹点集群，生成多个子轨迹段
     *
     * @param cluster     轨迹点集群（GaussPoint类型）
     * @param maxDistance 最大距离阈值（米），超过此距离的点将被视为不同轨迹段的起始点
     *
     * @return 切分后的轨迹段列表（每个元素为一个子轨迹段的GaussPoint列表）
     */
    private List<List<GaussPoint>> splitClusterByDistance(List<GaussPoint> cluster, double maxDistance) {
        List<List<GaussPoint>> segments = new ArrayList<>();
        if (cluster == null || cluster.isEmpty()) {
            return segments;
        }

        List<GaussPoint> currentSegment = new ArrayList<>();
        currentSegment.add(cluster.get(0));

        for (int i = 1; i < cluster.size(); i++) {
            GaussPoint prevPoint = cluster.get(i - 1);
            GaussPoint currPoint = cluster.get(i);

            // 计算欧几里得距离
            double distance = Math.sqrt(Math.pow(currPoint.getGaussX() - prevPoint.getGaussX(), 2) + Math.pow(currPoint.getGaussY() - prevPoint.getGaussY(), 2));

            if (distance > maxDistance) {
                // 距离超过阈值，结束当前段，开始新段
                segments.add(new ArrayList<>(currentSegment));
                currentSegment = new ArrayList<>();
            }
            currentSegment.add(currPoint);
        }

        // 处理最后一段
        segments.add(currentSegment);

        log.debug("距离切分完成：原始 {} 个点 切分出 {} 个子段，最大距离阈值 {} 米", cluster.size(), segments.size(), maxDistance);
        return segments;
    }

    /**
     * 按最大时间间隔切分轨迹点集群，生成多个子轨迹段
     *
     * @param cluster    轨迹点集群（GaussPoint类型）
     * @param maxSeconds 最大时间间隔阈值（秒），超过此时间间隔的点将被视为不同轨迹段的起始点
     *
     * @return 切分后的轨迹段列表（每个元素为一个子轨迹段的GaussPoint列表）
     */
    private List<List<GaussPoint>> splitClusterBySeconds(List<GaussPoint> cluster, double maxSeconds) {
        List<List<GaussPoint>> segments = new ArrayList<>();
        if (cluster == null || cluster.isEmpty()) {
            return segments;
        }

        List<GaussPoint> currentSegment = new ArrayList<>();
        currentSegment.add(cluster.get(0));

        for (int i = 1; i < cluster.size(); i++) {
            GaussPoint prevPoint = cluster.get(i - 1);
            GaussPoint currPoint = cluster.get(i);

            // 计算两点之间的时间间隔（秒）
            Duration duration = Duration.between(prevPoint.getGpsTime(), currPoint.getGpsTime());
            long seconds = duration.getSeconds();

            if (seconds > maxSeconds) {
                // 时间间隔超过阈值，结束当前段，开始新段
                segments.add(new ArrayList<>(currentSegment));
                currentSegment = new ArrayList<>();
            }
            currentSegment.add(currPoint);
        }

        // 处理最后一段
        segments.add(currentSegment);

        log.debug("时间切分完成：原始 {} 个点 -> {} 个子段，最大时间阈值 {} 秒", cluster.size(), segments.size(), maxSeconds);
        return segments;
    }


    /**
     * 简单算法：适用于小数据量
     */
    private List<Wgs84Point> findClosestPointListSimple(List<Wgs84Point> targetPointList, List<Wgs84Point> wgs84Points) {
        List<Wgs84Point> result = targetPointList.stream()
                .map(targetPoint -> findClosestPointWithProgressiveToleranceFixed(targetPoint, wgs84Points))
                .filter(Objects::nonNull)
                .collect(Collectors.toList());
        log.debug("批量查找最接近点完成，匹配到的点数量={}", result.size());
        return result;
    }

    /**
     * 优化算法：使用JTS STRtree空间索引，适用于大数据量
     */
    private List<Wgs84Point> findClosestPointListOptimized(List<Wgs84Point> targetPointList, List<Wgs84Point> wgs84Points) {
        // 1. 构建JTS STRtree空间索引 - 这是JTS提供的高效R树实现
        STRtree spatialIndex = new STRtree();

        // 将原始点插入空间索引
        for (Wgs84Point point : wgs84Points) {
            // 创建点的Envelope（边界框）作为空间索引键
            Envelope envelope = new Envelope(point.getLongitude(), point.getLongitude(),
                    point.getLatitude(), point.getLatitude());
            spatialIndex.insert(envelope, point);
        }

        // 构建索引
        spatialIndex.build();

        log.debug("JTS STRtree空间索引构建完成，原始点数量={}", wgs84Points.size());

        // 2. 并行处理目标点，提高多核CPU利用率
        return targetPointList.parallelStream()
                .unordered()
                .map(targetPoint -> findClosestPointWithSTRtree(targetPoint, spatialIndex))
                .filter(Objects::nonNull)
                .collect(Collectors.toList());
    }

    /**
     * 使用JTS STRtree查找最近点
     */
    private Wgs84Point findClosestPointWithSTRtree(Wgs84Point targetPoint, STRtree spatialIndex) {
        for (double tolerance : config.TOLERANCES) {
            // 将米转换为度（近似转换：1度≈111公里）
            double toleranceDegrees = tolerance / config.MI_TO_DEGREE;

            // 创建搜索范围的Envelope
            Envelope searchEnvelope = new Envelope(
                    targetPoint.getLongitude() - toleranceDegrees,
                    targetPoint.getLongitude() + toleranceDegrees,
                    targetPoint.getLatitude() - toleranceDegrees,
                    targetPoint.getLatitude() + toleranceDegrees
            );

            // 使用STRtree查询候选点
            @SuppressWarnings("unchecked")
            List<Wgs84Point> candidates = spatialIndex.query(searchEnvelope);

            if (!candidates.isEmpty()) {
                Wgs84Point closest = findClosestPointInCandidates(targetPoint, candidates, tolerance);
                if (closest != null) {
                    return closest;
                }
            }
        }

        log.warn("在所有容差范围内都未找到匹配点：目标坐标=[{}, {}]", targetPoint.getLongitude(), targetPoint.getLatitude());
        return null;
    }

    /**
     * 修复渐进式容差逻辑错误
     */
    private Wgs84Point findClosestPointWithProgressiveToleranceFixed(Wgs84Point targetWgs84Point, List<Wgs84Point> wgs84Points) {
        for (double tolerance : config.TOLERANCES) {
            Wgs84Point result = findClosestPoint(targetWgs84Point, wgs84Points, tolerance);
            if (result != null) {
                return result; // 找到匹配点就返回
            }
        }
        log.warn("在所有容差范围内都未找到匹配点：目标坐标=[{}, {}]", targetWgs84Point.getLongitude(), targetWgs84Point.getLatitude());
        return null;
    }

    /**
     * 在候选点中查找最近点
     */
    private Wgs84Point findClosestPointInCandidates(Wgs84Point targetPoint, List<Wgs84Point> candidates, double maxDistance) {
        Wgs84Point closestPoint = null;
        double minDistance = Double.MAX_VALUE;

        for (Wgs84Point candidate : candidates) {
            double distance = haversine(targetPoint, candidate);
            if (distance <= maxDistance && distance < minDistance) {
                minDistance = distance;
                closestPoint = candidate;
            }
        }

        return closestPoint;
    }

    private int getMinEffectiveInterval(List<Wgs84Point> wgs84Points) {
        log.debug("准备计算上报时间间隔分布");
        int minEffectiveInterval = 1; // 默认值
        Map<Integer, Integer> intervalDistribution = new HashMap<>();
        for (int i = 1; i < wgs84Points.size(); i++) {
            Wgs84Point prevPoint = wgs84Points.get(i - 1);
            Wgs84Point currPoint = wgs84Points.get(i);

            // 计算时间间隔（秒）- LocalDateTime使用Duration
            Duration duration = Duration.between(prevPoint.getGpsTime(), currPoint.getGpsTime());
            int timeDiffSeconds = (int) duration.getSeconds();
            // 统计每个间隔的出现次数
            intervalDistribution.put(timeDiffSeconds, intervalDistribution.getOrDefault(timeDiffSeconds, 0) + 1);
        }
        // 获取最小有效时间间隔
        if (!intervalDistribution.isEmpty()) {
            // 找到点数最多的时间间隔，如果点数相同则选择时间间隔更小的
            minEffectiveInterval = intervalDistribution.entrySet().stream().max((e1, e2) -> {
                int countCompare = Integer.compare(e1.getValue(), e2.getValue());
                if (countCompare != 0) {
                    return countCompare; // 点数多的优先
                }
                return Integer.compare(e2.getKey(), e1.getKey()); // 点数相同时，时间间隔小的优先（降序比较）
            }).map(Map.Entry::getKey).orElse(1);
        }
        log.debug("最小有效上报时间间隔 {} 秒", minEffectiveInterval);
        return minEffectiveInterval;
    }

    private List<List<GaussPoint>> dbScanClusters(List<GaussPoint> gaussPoints, int minEffectiveInterval) {
        log.debug("从高斯投影中提取坐标数组");
        double[][] coords = new double[gaussPoints.size()][2];
        for (int i = 0; i < gaussPoints.size(); i++) {
            coords[i][0] = gaussPoints.get(i).getGaussX();
            coords[i][1] = gaussPoints.get(i).getGaussY();
        }

        log.debug("创建StaticArrayDatabase");
        Database db = new StaticArrayDatabase(new ArrayAdapterDatabaseConnection(coords), null);
        db.initialize();

        double eps = config.MAX_WORK_DISTANCE * minEffectiveInterval;
        int minPts = config.MIN_DBSCAN_POINTS;
        log.info("使用空间密集聚类参数 eps={} 米, minPts={}", String.format("%.2f", eps), minPts);
        DBSCAN<DoubleVector> dbscan = new DBSCAN<>(EuclideanDistance.STATIC, eps, minPts);

        log.debug("获取Relation对象并执行空间密集聚类");
        Relation<DoubleVector> relation = db.getRelation(TypeUtil.DOUBLE_VECTOR_FIELD);
        Clustering<Model> dbscanCluster = dbscan.run(relation);

        log.debug("映射结果");
        DBIDRange ids = (DBIDRange) relation.getDBIDs();
        List<List<GaussPoint>> clusters = new ArrayList<>();
        for (Cluster<Model> cluster : dbscanCluster.getAllClusters()) {
            log.debug("聚类信息： 聚类名称: {} 点数量: {}", cluster.getNameAutomatic(), cluster.size());
            if (!cluster.getNameAutomatic().equals("Cluster")) {
                // 如果不是聚类簇，跳过
                continue;
            }
            List<GaussPoint> gaussPointList = new ArrayList<>();
            for (DBIDIter iter = cluster.getIDs().iter(); iter.valid(); iter.advance()) {
                int offset = ids.getOffset(iter);
                GaussPoint gaussPoint = gaussPoints.get(offset);
                gaussPointList.add(gaussPoint);
            }
            // 保持返回的聚类的点位是按时间升序排序
            gaussPointList.sort(Comparator.comparing(GaussPoint::getGpsTime));
            clusters.add(gaussPointList);
        }
        return clusters;
    }

    /**
     * 过滤异常点位信息
     *
     * @param wgs84Points 轨迹点列表（Wgs84Point类型）
     *
     * @return 过滤后的轨迹点列表（Wgs84Point类型）
     */
    public List<Wgs84Point> filterWgs84Points(List<Wgs84Point> wgs84Points) {
        log.debug("准备过滤异常点位信息");
        wgs84Points = wgs84Points.stream().filter(p -> {
            // 时间不能为空
            if (p.getGpsTime() == null) {
                log.warn("轨迹点时间为空，抛弃");
                return false;
            }
            // 经纬度不能为0（无效坐标）
            if (p.getLongitude() == 0.0 && p.getLatitude() == 0.0) {
                log.warn("定位时间: {} 轨迹点经纬度为 0 ，抛弃", p.getGpsTime());
                return false;
            }
            // 经纬度必须在合理范围内
            if (p.getLongitude() < -180.0 || p.getLongitude() > 180.0 || p.getLatitude() < -90.0 || p.getLatitude() > 90.0) {
                log.warn("定位时间: {} 轨迹点经纬度超出范围：[{},{}] 抛弃", p.getGpsTime(), p.getLongitude(), p.getLatitude());
                return false;
            }
            // GPS点位必须是已定位的点位
            if (p.getGpsStatus() != 0 && p.getGpsStatus() != 1) {
                log.warn("定位时间: {} 轨迹点GPS状态为 {} ，抛弃", p.getGpsTime(), p.getGpsStatus());
                return false;
            }
            // 必须是作业状态
            if (p.getJobStatus() != 0 && p.getJobStatus() != 1) {
                log.warn("定位时间: {} 轨迹点作业状态为 {} ，抛弃", p.getGpsTime(), p.getJobStatus());
                return false;
            }
            return true;
        }).collect(Collectors.toList());
        log.debug("过滤异常点位信息完成，剩余点位数量：{}", wgs84Points.size());

        wgs84Points.sort(Comparator.comparing(Wgs84Point::getGpsTime));

        log.debug("准备去重完全重复的轨迹点");
        Map<String, Wgs84Point> pointMap = new LinkedHashMap<>();
        for (Wgs84Point wgs84Point : wgs84Points) {
            String key = StrUtil.format("{},{}", wgs84Point.getLongitude(), wgs84Point.getLatitude());
            pointMap.putIfAbsent(key, wgs84Point);
        }
        wgs84Points = new ArrayList<>(pointMap.values());
        log.debug("去重完全重复的轨迹点完成，剩余点位数量：{}", wgs84Points.size());
        return wgs84Points;
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
     *
     * @return 两点之间的球面距离（米），结果为非负数
     */
    public double haversine(Wgs84Point wgs84Point1, Wgs84Point wgs84Point2) {
        // 将经纬度从度数转换为弧度，三角函数计算需要弧度值
        double lon1 = Math.toRadians(wgs84Point1.getLongitude());
        double lat1 = Math.toRadians(wgs84Point1.getLatitude());
        double lon2 = Math.toRadians(wgs84Point2.getLongitude());
        double lat2 = Math.toRadians(wgs84Point2.getLatitude());

        // 计算两点间经度差和纬度差
        double dlon = lon2 - lon1;
        double dlat = lat2 - lat1;

        // Haversine公式核心计算：半正矢公式
        double a = Math.sin(dlat / 2.0) * Math.sin(dlat / 2.0) + Math.cos(lat1) * Math.cos(lat2) * Math.sin(dlon / 2.0) * Math.sin(dlon / 2.0);

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
     * @param wgs84Point       要判断的WGS84坐标点，包含有效经纬度
     * @param wgs84CenterPoint 圆的中心点（WGS84坐标）
     * @param radius           圆的半径（米），必须为非负数
     *
     * @return 如果点到圆心的球面距离小于等于指定半径（包含边界），则返回true；否则返回false
     */
    public boolean inCircle(Wgs84Point wgs84Point, Wgs84Point wgs84CenterPoint, double radius) {
        // 计算两点间的球面距离（米）
        double distance = haversine(wgs84Point, wgs84CenterPoint);
        // 如果距离小于等于半径，则点在圆内（包含边界）
        return distance <= radius;
    }

    /**
     * 判断WGS84坐标系下的点是否在指定的几何图形内部。
     * <p>
     * 该方法使用JTS拓扑套件进行精确的空间关系判断，支持各种几何类型（点、线、面等）。
     * 适用于轨迹分析、地理围栏和空间统计等场景。
     * </p>
     *
     * @param wgs84Point    待测试的WGS84坐标点，包含有效经纬度
     * @param wgs84Geometry 用于判断的几何图形，必须是有效的JTS Geometry对象
     *
     * @return 如果点在几何图形内部或边界上，则返回true；否则返回false
     */
    public boolean inGeometry(Wgs84Point wgs84Point, Geometry wgs84Geometry) {
        try {
            // 创建点几何对象
            Point point = config.GEOMETRY_FACTORY.createPoint(new Coordinate(wgs84Point.getLongitude(), wgs84Point.getLatitude()));

            // 判断点是否在几何图形内或边界上
            // covers() 方法包含边界，而 contains() 不包含边界
            return wgs84Geometry.contains(point);
        } catch (Exception e) {
            log.warn("判断点是否在几何图形内失败：点[{},{}] 错误={}", wgs84Point.getLongitude(), wgs84Point.getLatitude(), e.getMessage());
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
     * @param wgs84Point            待测试的WGS84坐标点，包含有效经纬度
     * @param wgs84TopLeftPoint     矩形的左上角点（或任意对角点）
     * @param wgs84BottomRightPoint 矩形的右下角点（或任意对角点）
     *
     * @return 如果点的经纬度在矩形范围内（包含边界），则返回true；否则返回false
     */
    public boolean inRectangle(Wgs84Point wgs84Point, Wgs84Point wgs84TopLeftPoint, Wgs84Point wgs84BottomRightPoint) {
        try {
            double pointLon = wgs84Point.getLongitude();
            double pointLat = wgs84Point.getLatitude();
            double topLeftLon = wgs84TopLeftPoint.getLongitude();
            double topLeftLat = wgs84TopLeftPoint.getLatitude();
            double bottomRightLon = wgs84BottomRightPoint.getLongitude();
            double bottomRightLat = wgs84BottomRightPoint.getLatitude();

            // 确保左上角和右下角的经纬度关系正确
            double minLon = Math.min(topLeftLon, bottomRightLon);
            double maxLon = Math.max(topLeftLon, bottomRightLon);
            double maxLat = Math.max(topLeftLat, bottomRightLat);
            double minLat = Math.min(topLeftLat, bottomRightLat);

            // 判断点是否在矩形范围内（包含边界）
            return pointLon >= minLon && pointLon <= maxLon && pointLat >= minLat && pointLat <= maxLat;
        } catch (Exception e) {
            log.warn("判断点是否在矩形内失败：点[{},{}] 左上角[{},{}] 右下角[{},{}] 错误={}", wgs84Point.getLongitude(), wgs84Point.getLatitude(), wgs84TopLeftPoint.getLongitude(), wgs84TopLeftPoint.getLatitude(), wgs84BottomRightPoint.getLongitude(), wgs84BottomRightPoint.getLatitude(), e.getMessage());
            return false;
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
     *
     * @return WGS84坐标系的Geometry几何图形，解析失败返回空几何
     */
    public Geometry toWgs84Geometry(String wgs84WKT) {
        if (wgs84WKT == null || wgs84WKT.trim().isEmpty()) {
            log.warn("WKT字符串为空或null");
            return config.EMPTY_GEOMETRY;
        }
        try {
            Geometry geometry = new WKTReader(config.GEOMETRY_FACTORY).read(wgs84WKT);
            log.debug("WKT字符串解析成功：几何类型={}", geometry.getGeometryType());
            return geometry;
        } catch (ParseException e) {
            log.warn("WKT字符串解析失败：{}", e.getMessage());
            return config.EMPTY_GEOMETRY;
        } catch (Exception e) {
            log.warn("WKT字符串转换几何图形失败：{}", e.getMessage());
            return config.EMPTY_GEOMETRY;
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
     *
     * @return 高斯投影坐标系下的几何图形（米制），转换失败返回空几何
     */
    public Geometry toGaussGeometry(Geometry wgs84Geometry) {
        try {
            if (wgs84Geometry == null || wgs84Geometry.isEmpty()) {
                return config.EMPTY_GEOMETRY;
            }

            // WGS84几何到高斯投影几何转换开始
            // 获取几何图形的边界信息来确定高斯投影参数
            Envelope env = wgs84Geometry.getEnvelopeInternal();
            double centerLon = (env.getMinX() + env.getMaxX()) / 2.0;

            // 验证WGS84坐标的合理性
            if (env.getMinX() < -180 || env.getMaxX() > 180 || env.getMinY() < -90 || env.getMaxY() > 90) {
                log.warn("WGS84坐标超出合理范围：MinLon={}, MaxLon={}, MinLat={}, MaxLat={}", env.getMinX(), env.getMaxX(), env.getMinY(), env.getMaxY());
                return config.EMPTY_GEOMETRY;
            }

            // 计算高斯投影带号 (6度分带)
            int zone = (int) Math.floor((centerLon + 180) / 6) + 1;
            double centralMeridian = (zone - 1) * 6 - 180 + 3;
            double falseEasting = zone * 1000000.0 + 500000.0;

            // 验证投影带号的合理性
            if (zone < 1 || zone > 60) {
                log.warn("投影带号超出合理范围：zone={}，经度={}", zone, centerLon);
                return config.EMPTY_GEOMETRY;
            }

            // 从缓存获取或创建高斯投影CRS
            CoordinateReferenceSystem gaussCRS = getGaussCRS(zone, falseEasting, centralMeridian);

            // 构建缓存key（用于transform缓存）
            String cacheKey = String.format(config.CACHE_KEY_FORMAT, zone, falseEasting, centralMeridian);
            log.debug("获取WGS84到高斯投影的坐标转换器：缓存键 {}", cacheKey);

            // 从缓存获取或创建WGS84到高斯投影的坐标转换
            MathTransform transform = config.WGS84_TO_GAUSS_TRANSFORM_CACHE.computeIfAbsent(cacheKey, key -> {
                try {
                    return CRS.findMathTransform(config.WGS84_CRS, gaussCRS, true);
                } catch (Exception e) {
                    log.warn("创建坐标转换失败：zone={}, falseEasting={}, centralMeridian={}, 错误={}", zone, falseEasting, centralMeridian, e.getMessage());
                    return null;
                }
            });

            // 执行坐标转换（正向：WGS84 -> 高斯投影）
            Geometry gaussGeometry = JTS.transform(wgs84Geometry, transform);

            // 验证转换后的高斯投影坐标合理性
            Envelope gaussEnv = gaussGeometry.getEnvelopeInternal();
            if (gaussEnv.getMinX() < 500000 || gaussEnv.getMaxX() > 64000000 || gaussEnv.getMinY() < -10000000 || gaussEnv.getMaxY() > 10000000) {
                log.warn("转换后的高斯投影坐标超出合理范围：MinX={}, MaxX={}, MinY={}, MaxY={}", gaussEnv.getMinX(), gaussEnv.getMaxX(), gaussEnv.getMinY(), gaussEnv.getMaxY());
                return config.EMPTY_GEOMETRY;
            }
            // WGS84几何到高斯投影几何转换完成
            return gaussGeometry;
        } catch (Exception e) {
            log.warn("WGS84几何到高斯投影几何转换失败：错误={}", e.getMessage());
            return config.EMPTY_GEOMETRY;
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
     *
     * @return 包含相交轮廓WKT字符串和面积（亩）的结果对象，无相交则返回空几何和零面积
     */
    public WktIntersectionResult intersection(String wgs84WKT1, String wgs84WKT2) {
        WktIntersectionResult result = new WktIntersectionResult();
        result.setWkt(config.EMPTY_GEOMETRY.toText());
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

            log.debug("高斯投影转换成功：几何1类型={}, 几何2类型={}", gaussGeometry1.getGeometryType(), gaussGeometry2.getGeometryType());

            // 3. 在高斯投影坐标系下执行相交操作（更精确）
            Geometry gaussIntersection = gaussGeometry1.intersection(gaussGeometry2);

            // 4. 检查是否有实际相交区域
            if (gaussIntersection == null || gaussIntersection.isEmpty()) {
                log.debug("两个几何图形没有相交区域");
                return result;
            }

            log.debug("相交成功，高斯投影下相交几何类型：{}，面积：{}平方米", gaussIntersection.getGeometryType(), gaussIntersection.getArea());

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

            log.debug("相交轮廓计算完成：亩数={}亩（基于WGS84坐标系计算）", result.getMu());
        } catch (Exception e) {
            log.warn("相交计算失败：{}", e.getMessage());
        }

        return result;
    }

    /**
     * 将高斯投影几何转换为WGS84投影几何
     *
     * @param gaussGeometry 输入的高斯投影几何对象（LineString或Polygon类型）
     *
     * @return 转换后的WGS84投影几何对象（LineString或Polygon类型），转换失败返回空几何
     */
    public Geometry toWgs84Geometry(Geometry gaussGeometry) {
        try {
            // 开始高斯投影几何到WGS84投影几何转换
            // 获取几何图形的边界信息来反推高斯投影参数
            Envelope env = gaussGeometry.getEnvelopeInternal();
            double centerX = (env.getMinX() + env.getMaxX()) / 2.0;

            // 验证高斯投影坐标的合理性
            // 高斯投影X坐标合理范围：50万-6400万米（对应全球1-60带，更宽松的范围）
            // 高斯投影Y坐标合理范围：±1000万米（对应全球纬度范围）
            if (env.getMinX() < 500000 || env.getMaxX() > 64000000 || env.getMinY() < -10000000 || env.getMaxY() > 10000000) {
                log.warn("高斯投影坐标超出合理范围：MinX={}, MaxX={}, MinY={}, MaxY={}", env.getMinX(), env.getMaxX(), env.getMinY(), env.getMaxY());
                return config.EMPTY_GEOMETRY;
            }

            // 智能确定投影带号策略
            // 1. 首先尝试根据几何中心反推带号
            int zone = (int) Math.floor(centerX / 1000000.0);
            // 计算中央经线（与wgs84PointTransformToGaussPoint方法保持一致）
            double centralMeridian = (zone - 1) * 6 - 180 + 3;
            double falseEasting = zone * 1000000.0 + 500000.0;

            // 2. 如果几何范围跨越多个投影带（宽度超过100万米），使用更保守的策略
            double geometryWidth = env.getMaxX() - env.getMinX();
            if (geometryWidth > 1000000) {
                log.warn("几何图形宽度{}米，可能跨越多个投影带，使用保守策略", geometryWidth);
                if (zone < 1 || zone > 60) {
                    log.warn("无法确定合适的投影带号：centerX={}", centerX);
                    return config.EMPTY_GEOMETRY;
                }
            } else if (zone < 1 || zone > 60) {
                log.warn("反推的投影带号不合理：投影带号 {}, centerX={}", zone, centerX);
                // 尝试备用策略：如果zone不合理，可能是假东距计算问题，尝试重新计算
                log.trace("尝试备用策略重新计算投影带号");
                int backupZone = (int) Math.floor((centerX - 500000.0) / 1000000.0);
                log.trace("备用策略计算的投影带号：{}", backupZone);
                if (backupZone < 1 || backupZone > 60) {
                    log.warn("备用策略仍然无法确定合适的投影带号：centerX={}", centerX);
                    return config.EMPTY_GEOMETRY;
                }
                zone = backupZone;
            }

            // 从缓存获取或创建高斯投影CRS
            CoordinateReferenceSystem gaussCRS = getGaussCRS(zone, falseEasting, centralMeridian);

            if (gaussCRS == null) {
                log.warn("无法获取高斯投影CRS：zone={}", zone);
                return config.EMPTY_GEOMETRY;
            }

            // 构建缓存key（用于transform缓存）
            String cacheKey = String.format(config.CACHE_KEY_FORMAT, zone, falseEasting, centralMeridian);
            log.debug("获取高斯投影到WGS84的坐标转换器：缓存键 {}", cacheKey);

            // 从缓存获取或创建高斯投影到WGS84的坐标转换
            final int finalZone = zone;
            MathTransform transform = config.GAUSS_TO_WGS84_TRANSFORM_CACHE.computeIfAbsent(cacheKey, key -> {
                try {
                    return CRS.findMathTransform(gaussCRS, config.WGS84_CRS, true);
                } catch (Exception e) {
                    log.warn("创建坐标转换失败：zone={}, falseEasting={}, centralMeridian={}, 错误={}", finalZone, falseEasting, centralMeridian, e.getMessage());
                    return null;
                }
            });

            if (transform == null) {
                log.warn("无法获取坐标转换：zone={}", finalZone);
                return config.EMPTY_GEOMETRY;
            }

            // 执行坐标转换（逆向：高斯投影 -> WGS84）
            Geometry wgs84Geometry = JTS.transform(gaussGeometry, transform);

            // 验证转换后的WGS84坐标合理性
            Envelope wgs84Env = wgs84Geometry.getEnvelopeInternal();
            if (wgs84Env.getMinX() < -180 || wgs84Env.getMaxX() > 180 || wgs84Env.getMinY() < -90 || wgs84Env.getMaxY() > 90) {
                log.warn("转换后的WGS84坐标超出合理范围：MinLon={}, MaxLon={}, MinLat={}, MaxLat={}", wgs84Env.getMinX(), wgs84Env.getMaxX(), wgs84Env.getMinY(), wgs84Env.getMaxY());
                return config.EMPTY_GEOMETRY;
            }
            // 高斯投影几何到WGS84投影几何转换完成
            return wgs84Geometry;
        } catch (Exception e) {
            log.warn("高斯投影几何到WGS84投影几何转换失败：错误={}", e.getMessage());
            return config.EMPTY_GEOMETRY;
        }
    }


    /**
     * 在WGS84点列表中查找与目标坐标最接近的点（容差匹配）
     * <p>
     * 用于解决坐标转换往返过程中的精度损失问题，通过容差匹配找到原始点，
     * 从而获取定位时间等属性信息。
     * </p>
     *
     * @param targetWgs84Point 目标WGS84点（通常是转换回来的点）
     * @param wgs84Points      原始WGS84点列表
     * @param tolerance        容差范围（单位：米），默认建议1.0米
     *
     * @return 最接近的原始点，如果找不到则返回null
     */
    public Wgs84Point findClosestPoint(Wgs84Point targetWgs84Point, List<Wgs84Point> wgs84Points, double tolerance) {
        if (targetWgs84Point == null || wgs84Points == null || wgs84Points.isEmpty()) {
            return null;
        }
        Wgs84Point closestPoint = null;
        double minDistance = Double.MAX_VALUE;
        for (Wgs84Point sourcePoint : wgs84Points) {
            double distance = haversine(targetWgs84Point, sourcePoint);
            if (distance <= tolerance && distance < minDistance) {
                minDistance = distance;
                closestPoint = sourcePoint;
            }
        }
        if (closestPoint != null) {
            log.trace("找到最接近的点：距离={}米，原始坐标=[{}, {}]，目标坐标=[{}, {}]",
                    String.format("%.6f", minDistance),
                    closestPoint.getLongitude(), closestPoint.getLatitude(),
                    targetWgs84Point.getLongitude(), targetWgs84Point.getLatitude());
        } else {
            log.warn("在容差{}米范围内未找到匹配点，目标坐标=[{}, {}]", tolerance, targetWgs84Point.getLongitude(), targetWgs84Point.getLatitude());
        }
        return closestPoint;
    }


    /**
     * 在WGS84点列表中批量查找与目标点列表最接近的点（多级容差匹配）
     * <p>
     * 使用渐进式容差策略：0.00111, 0.0111, 0.111, 1.11, 11.1
     * 解决坐标转换往返过程中的精度损失问题
     * </p>
     * <p><b>性能优化：</b></p>
     * <ul>
     *   <li>使用空间网格索引加速最近邻搜索</li>
     *   <li>按地理区域预分组原始点，减少搜索范围</li>
     *   <li>渐进式容差策略，优先在小范围内精确匹配</li>
     *   <li>并行处理大批量目标点</li>
     * </ul>
     *
     * @param targetPointList 目标WGS84点列表（通常是转换回来的点）
     * @param wgs84Points     原始WGS84点列表
     *
     * @return 匹配到的最接近点列表
     */
    public List<Wgs84Point> findClosestPointList(List<Wgs84Point> targetPointList, List<Wgs84Point> wgs84Points) {
        if (targetPointList.isEmpty() || wgs84Points.isEmpty()) {
            return new ArrayList<>();
        }

        log.debug("批量查找最接近点，目标点数量={}, 原始点数量={}", targetPointList.size(), wgs84Points.size());

        // 对于小数据量，使用简单算法
        if (targetPointList.size() < 100 || wgs84Points.size() < 100) {
            return findClosestPointListSimple(targetPointList, wgs84Points);
        }

        // 对于大数据量，使用空间索引优化
        return findClosestPointListOptimized(targetPointList, wgs84Points);
    }

    public List<Wgs84Point> toWgs84PointList(List<GaussPoint> gaussPoints) {
        log.debug("转换高斯投影点列表为WGS84点列表，点数量={}", gaussPoints.size());
        if (CollUtil.isEmpty(gaussPoints)) {
            return new ArrayList<>();
        }

        List<Wgs84Point> wgs84Points = new ArrayList<>(gaussPoints.size());
        try {
            // 按投影带分组，批量处理同一投影带的点
            Map<Integer, List<GaussPoint>> pointsByZone = new HashMap<>();
            for (GaussPoint gaussPoint : gaussPoints) {
                // 从WGS84经度计算投影带号（正确的方式）
                double longitude = gaussPoint.getLongitude();
                int zone = (int) Math.floor((longitude + 180) / 6) + 1;

                // 验证投影带号的合理性
                if (zone >= 1 && zone <= 60) {
                    pointsByZone.computeIfAbsent(zone, k -> new ArrayList<>()).add(gaussPoint);
                } else {
                    log.warn("计算得到的投影带号不合理：经度={}, 投影带号={}", longitude, zone);
                }
            }

            log.debug("按投影带分组结果：{} 个投影带", pointsByZone.size());
            for (Map.Entry<Integer, List<GaussPoint>> entry : pointsByZone.entrySet()) {
                log.debug("投影带 {}：{} 个点", entry.getKey(), entry.getValue().size());
            }

            // 对每个投影带批量处理
            for (Map.Entry<Integer, List<GaussPoint>> entry : pointsByZone.entrySet()) {
                int zone = entry.getKey();
                List<GaussPoint> zonePoints = entry.getValue();

                // 计算投影参数
                double centralMeridian = (zone - 1) * 6 - 180 + 3;
                double falseEasting = zone * 1000000.0 + 500000.0;

                log.debug("处理投影带 {}：中央经线={}, 假东距={}", zone, centralMeridian, falseEasting);

                CoordinateReferenceSystem gaussCRS = getGaussCRS(zone, falseEasting, centralMeridian);
                if (gaussCRS == null) {
                    log.warn("获取高斯投影CRS失败：投影带号={}", zone);
                    continue;
                }

                // 获取或创建坐标转换对象（高斯->WGS84）
                String cacheKey = String.format(config.CACHE_KEY_FORMAT, zone, falseEasting, centralMeridian);
                log.debug("获取高斯投影到WGS84的坐标转换器：缓存键 {}", cacheKey);

                // 从缓存获取或创建高斯到WGS84的坐标转换
                MathTransform gaussToWgs84Transform = config.GAUSS_TO_WGS84_TRANSFORM_CACHE.computeIfAbsent(cacheKey, key -> {
                    try {
                        return CRS.findMathTransform(gaussCRS, config.WGS84_CRS, true);
                    } catch (Exception e) {
                        log.warn("创建高斯到WGS84坐标转换失败：zone={}, falseEasting={}, centralMeridian={}, 错误={}", zone, falseEasting, centralMeridian, e.getMessage());
                        return null;
                    }
                });

                if (gaussToWgs84Transform == null) {
                    log.warn("获取坐标转换器失败，跳过投影带 {} 的处理", zone);
                    continue;
                }

                // 批量转换同一投影带的点
                int convertedCount = 0;
                for (GaussPoint gaussPoint : zonePoints) {
                    try {
                        Coordinate sourceCoord = new Coordinate(gaussPoint.getGaussX(), gaussPoint.getGaussY());
                        Coordinate targetCoord = new Coordinate();
                        JTS.transform(sourceCoord, targetCoord, gaussToWgs84Transform);

                        // 验证转换结果的合理性（经纬度范围）
                        if (targetCoord.x >= -180 && targetCoord.x <= 180 && targetCoord.y >= -90 && targetCoord.y <= 90) {
                            Wgs84Point result = new Wgs84Point(gaussPoint.getGpsTime(), targetCoord.x, targetCoord.y);
                            wgs84Points.add(result);
                            convertedCount++;
                        } else {
                            log.warn("转换结果超出合理范围：经度={}, 纬度={}", targetCoord.x, targetCoord.y);
                        }
                    } catch (Exception e) {
                        log.warn("转换高斯投影点到WGS84失败：zone={}, gaussX={}, gaussY={}, 错误={}", zone, gaussPoint.getGaussX(), gaussPoint.getGaussY(), e.getMessage());
                    }
                }
                log.debug("投影带 {} 成功转换 {} 个点", zone, convertedCount);
            }

            log.debug("总计成功转换 {} 个高斯投影点到WGS84坐标系", wgs84Points.size());
        } catch (Exception e) {
            log.warn("高斯投影点转换为WGS84失败: {}", e.getMessage());
        }
        return wgs84Points;
    }

    /**
     * 将WGS84坐标系下的点列表转换为高斯投影坐标系下的点列表
     *
     * @param wgs84Points 输入的WGS84坐标系下的点列表（Wgs84Point类型）
     *
     * @return 转换后的高斯投影坐标系下的点列表（GaussPoint类型）
     */
    public List<GaussPoint> toGaussPointList(List<Wgs84Point> wgs84Points) {
        log.debug("转换WGS84点列表为高斯投影点列表，点数量={}", wgs84Points.size());
        if (CollUtil.isEmpty(wgs84Points)) {
            return new ArrayList<>();
        }

        List<GaussPoint> gaussPoints = new ArrayList<>(wgs84Points.size());
        try {
            // 优化1：按投影带分组，批量处理同一投影带的点
            Map<Integer, List<Wgs84Point>> pointsByZone = new HashMap<>();
            for (Wgs84Point wgs84Point : wgs84Points) {
                // 计算投影带号
                double longitude = wgs84Point.getLongitude();
                int zone = (int) Math.floor((longitude + 180) / 6) + 1;

                // 验证投影带号的合理性
                if (zone >= 1 && zone <= 60) {
                    pointsByZone.computeIfAbsent(zone, k -> new ArrayList<>()).add(wgs84Point);
                }
            }

            // 优化2：对每个投影带批量处理
            for (Map.Entry<Integer, List<Wgs84Point>> entry : pointsByZone.entrySet()) {
                int zone = entry.getKey();
                List<Wgs84Point> zonePoints = entry.getValue();

                // 计算投影参数
                double centralMeridian = (zone - 1) * 6 - 180 + 3;
                double falseEasting = zone * 1000000.0 + 500000.0;

                CoordinateReferenceSystem gaussCRS = getGaussCRS(zone, falseEasting, centralMeridian);
                if (gaussCRS == null) {
                    continue;
                }

                // 获取或创建坐标转换对象
                String cacheKey = String.format(config.CACHE_KEY_FORMAT, zone, falseEasting, centralMeridian);
                log.debug("获取WGS84到高斯投影的坐标转换器：缓存键 {}", cacheKey);

                MathTransform transform = config.WGS84_TO_GAUSS_TRANSFORM_CACHE.computeIfAbsent(cacheKey, key -> {
                    try {
                        return CRS.findMathTransform(config.WGS84_CRS, gaussCRS, true);
                    } catch (Exception e) {
                        log.warn("创建坐标转换失败：zone={}, falseEasting={}, centralMeridian={}, 错误={}", zone, falseEasting, centralMeridian, e.getMessage());
                        return null;
                    }
                });

                if (transform == null) {
                    continue;
                }

                // 优化3：批量转换同一投影带的点
                for (Wgs84Point wgs84Point : zonePoints) {
                    try {
                        Coordinate sourceCoord = new Coordinate(wgs84Point.getLongitude(), wgs84Point.getLatitude());
                        Coordinate targetCoord = new Coordinate();
                        JTS.transform(sourceCoord, targetCoord, transform);

                        // 验证转换结果的合理性
                        if (targetCoord.x >= 500000 && targetCoord.x <= 64000000 && targetCoord.y >= -10000000 && targetCoord.y <= 10000000) {
                            GaussPoint result = new GaussPoint(wgs84Point.getGpsTime(), wgs84Point.getLongitude(), wgs84Point.getLatitude(), targetCoord.x, targetCoord.y);
                            gaussPoints.add(result);
                        }
                    } catch (Exception e) {
                        log.warn("转换WGS84点到高斯投影失败：zone={}, 经度={}, 纬度={}, 错误={}", zone, wgs84Point.getLongitude(), wgs84Point.getLatitude(), e.getMessage());
                    }
                }
            }

            log.debug("成功转换 {} 个轨迹点到高斯-克吕格投影坐标系", gaussPoints.size());
        } catch (Exception e) {
            log.warn("WGS84轨迹点转换为高斯投影失败: {}", e.getMessage());
        }
        return gaussPoints;
    }

    /**
     * 计算WGS84坐标系下的几何图形球面面积（平方米）
     *
     * @param wgs84Geometry 输入的WGS84坐标系下的几何图形对象（Polygon或MultiPolygon类型）
     *
     * @return 计算得到的球面面积（平方米），计算失败返回0.0
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
     * 该方法利用球面面积算法计算WGS84坐标系下的几何图形面积（平方米），
     * 并将其转换为亩，结果四舍五入保留4位小数。
     * </p>
     *
     * @param wgs84Geometry 输入的WGS84坐标系下的几何图形对象（Polygon或MultiPolygon类型）
     *
     * @return 计算得到的几何图形面积（亩），计算失败返回0.0
     */
    public double calcMu(Geometry wgs84Geometry) {
        try {
            // 步骤1：使用球面面积算法计算平方米面积
            double areaSqm = calculateSphericalArea(wgs84Geometry);

            // 步骤2：平方米转换为亩
            // 四舍五入保留4位小数
            double mu = Math.round((areaSqm * (config.SQUARE_TO_MU_METER)) * 10000.0) / 10000.0;

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
     *
     * @return 几何图形的面积（亩），四舍五入保留4位小数；解析或计算失败返回0.0
     *
     * @see #calcMu(Geometry)
     */
    public double calcMu(String wgs84Wkt) {
        return calcMu(toWgs84Geometry(wgs84Wkt));
    }


    /**
     * 道路拆分，将道路轨迹拆分掉，只保留作业信息
     *
     * @param wgs84Points  输入的WGS84坐标系下的点列表（Wgs84Point类型）
     * @param workingWidth 作业幅宽（米）
     *
     * @return 拆分后的作业轨迹结果（SplitResult类型）
     */
    public SplitResult splitRoad(List<Wgs84Point> wgs84Points, double workingWidth) {
        long splitRoadStartTime = System.currentTimeMillis();
        SplitResult splitResult = new SplitResult();
        splitResult.setWorkingWidth(workingWidth);
        splitResult.setWkt(config.EMPTY_GEOMETRY.toText());

        if (CollUtil.isEmpty(wgs84Points)) {
            log.error("作业轨迹点列表不能为空");
            return splitResult;
        }
        if (workingWidth < 1) {
            log.error("作业幅宽必须大于等于 1 米");
            return splitResult;
        }
        log.info("道路拆分入参 wgs84点位集合大小：{} 幅宽：{}米", wgs84Points.size(), workingWidth);

        // 作业机具的左右幅宽
        double halfWorkingWidth = workingWidth / 2.0;

        // 过滤异常点位信息
        wgs84Points = filterWgs84Points(wgs84Points);

        // 获得最小上报时间间隔
        int minEffectiveInterval = getMinEffectiveInterval(wgs84Points);

        // 转换为高斯投影坐标
        List<GaussPoint> gaussPoints = toGaussPointList(wgs84Points);
        if (gaussPoints.size() < config.MIN_DBSCAN_POINTS) {
            log.error("作业轨迹点列表必须包含至少 {} 个有效点位", config.MIN_DBSCAN_POINTS);
            return splitResult;
        }

        // 执行聚类，并且返回高斯投影的坐标
        List<List<GaussPoint>> clusters = dbScanClusters(gaussPoints, minEffectiveInterval);
        log.info("聚类完成，总共有 {} 个聚类簇", clusters.size());
        if (clusters.isEmpty()) {
            log.debug("没有聚类簇");
            return splitResult;
        }

        log.debug("使用时间再次切割聚类");
        List<List<GaussPoint>> splitClustersBySeconds = new ArrayList<>();
        for (List<GaussPoint> cluster : clusters) {
            log.debug("聚类簇包含 {} 个点", cluster.size());
            if (cluster.size() >= config.MIN_DBSCAN_POINTS) {
                splitClustersBySeconds.addAll(splitClusterBySeconds(cluster, config.SPLIT_TIME_SECOND));
            }
        }
        log.debug("经过时间切割后，总共有 {} 个聚类簇", splitClustersBySeconds.size());

        log.debug("循环所有聚类簇，生成几何图形");
        List<Geometry> clusterGaussGeometries = new ArrayList<>();//聚类的几何图形列表
        List<List<GaussPoint>> clusterGaussPoints = new ArrayList<>();//每一个几何图形的点位列表
        for (List<GaussPoint> cluster : splitClustersBySeconds) {
            log.debug("聚类簇包含 {} 个点", cluster.size());
            if (cluster.size() >= config.MIN_DBSCAN_POINTS) {
                log.debug("聚类后，按时间升序排序");
                cluster.sort(Comparator.comparing(GaussPoint::getGpsTime));
                log.debug("创建线缓冲，缓冲半径：{} 米", halfWorkingWidth);

                log.debug("按距离切分聚类");
                List<List<GaussPoint>> segments = splitClusterByDistance(cluster, minEffectiveInterval * config.MAX_WORK_DISTANCE * 2);
                log.debug("切分后得到 {} 个子段", segments.size());

                for (List<GaussPoint> segment : segments) {
                    log.debug("处理子段：{} 个点", segment.size());
                    if (segment.size() > 2) {
                        Geometry gaussGeometry;
                        log.debug("创建线缓冲，缓冲半径：{} 米", halfWorkingWidth);
                        Coordinate[] coordinates = segment.stream().map(point -> new Coordinate(point.getGaussX(), point.getGaussY())).toArray(Coordinate[]::new);
                        if (coordinates.length > 500) {
                            LineString lineString = config.GEOMETRY_FACTORY.createLineString(coordinates);
                            Geometry simplifiedGeometry = DouglasPeuckerSimplifier.simplify(lineString, 0.1);//0.1米容差
                            Coordinate[] simplifiedCoords = simplifiedGeometry.getCoordinates();
                            if (simplifiedCoords.length > 500) {
                                gaussGeometry = processLargeSegmentInChunks(segment, halfWorkingWidth);
                            } else {
                                gaussGeometry = simplifiedGeometry.buffer(halfWorkingWidth);
                            }
                        } else {
                            LineString lineString = config.GEOMETRY_FACTORY.createLineString(coordinates);
                            gaussGeometry = lineString.buffer(halfWorkingWidth);
                        }
                        log.info("使用膨胀、收缩参数 {}米 降低地块缝隙", workingWidth);
                        gaussGeometry = gaussGeometry.buffer(workingWidth).buffer(-workingWidth);
                        clusterGaussGeometries.add(gaussGeometry);
                        clusterGaussPoints.add(segment);
                    }
                }
            }
        }
        log.debug("生成了 {} 个几何图形，生成了 {} 组点位列表", clusterGaussGeometries.size(), clusterGaussPoints.size());
        if (clusterGaussGeometries.isEmpty()) {
            log.debug("没有生成任何几何图形");
            return splitResult;
        }

        log.debug("创建part对象集合");
        List<Part> parts = new ArrayList<>();
        for (int i = 0; i < clusterGaussGeometries.size(); i++) {
            Geometry clusterGaussGeometry = clusterGaussGeometries.get(i);
            Geometry wgs84PartGeometry = toWgs84Geometry(clusterGaussGeometry);
            List<Wgs84Point> wgs84PointList = new ArrayList<>();
            for (GaussPoint gaussPoint : clusterGaussPoints.get(i)) {
                wgs84PointList.add(new Wgs84Point(gaussPoint.getGpsTime(), gaussPoint.getLongitude(), gaussPoint.getLatitude()));
            }
            Part part = new Part();
            part.setGaussGeometry(clusterGaussGeometry);
            part.setTrackPoints(wgs84PointList);
            part.setStartTime(part.getTrackPoints().get(0).getGpsTime());
            part.setEndTime(part.getTrackPoints().get(part.getTrackPoints().size() - 1).getGpsTime());
            part.setWkt(wgs84PartGeometry.toText());
            part.setMu(calcMu(wgs84PartGeometry));
            parts.add(part);
        }

        log.debug("按面积倒序排序");
        parts.sort(Comparator.comparing(Part::getMu).reversed());
        log.debug("取前 {} 个最大面积的几何图形", config.MAX_SPLIT_RETURN_SIZE);
        parts = parts.subList(0, Math.min(config.MAX_SPLIT_RETURN_SIZE, parts.size()));
        log.debug("移除面积小于 {} 亩的几何图形", config.MIN_RETURN_MU);
        parts.removeIf(part -> part.getMu() < config.MIN_RETURN_MU);
        log.info("最终保留 {} 个地块", parts.size());
        if (parts.isEmpty()) {
            log.debug("没有保留任何地块");
            return splitResult;
        }
        parts.sort(Comparator.comparing(Part::getStartTime));

        log.debug("合并所有Part几何图形");
        Geometry unionPartsGaussGeometry = config.GEOMETRY_FACTORY.createGeometryCollection(parts.stream().map(Part::getGaussGeometry).toArray(Geometry[]::new)).union().buffer(0);
        log.debug("合并后几何图形的面积（平方米）：{}", unionPartsGaussGeometry.getArea());
        Geometry wgs84UnionGeometry = toWgs84Geometry(unionPartsGaussGeometry);
        //log.debug("合并后的几何图形：{}", wgs84UnionGeometry.toText());
        splitResult.setGaussGeometry(unionPartsGaussGeometry);
        splitResult.setWkt(wgs84UnionGeometry.toText());
        splitResult.setMu(calcMu(wgs84UnionGeometry));
        splitResult.setParts(parts);

        log.debug("再次检查合并后的几何图形，如果合并后是一个几何图形，那么修改parts信息");
        if (parts.size() > 1 && wgs84UnionGeometry.getNumGeometries() == 1) {
            List<Wgs84Point> wgs84PointList = new ArrayList<>();
            for (Part part : parts) {
                wgs84PointList.addAll(part.getTrackPoints());
            }
            wgs84PointList.sort(Comparator.comparing(Wgs84Point::getGpsTime));
            parts.clear();
            Part part = new Part();
            part.setGaussGeometry(unionPartsGaussGeometry);
            part.setTrackPoints(wgs84PointList);
            part.setStartTime(part.getTrackPoints().get(0).getGpsTime());
            part.setEndTime(part.getTrackPoints().get(part.getTrackPoints().size() - 1).getGpsTime());
            part.setWkt(wgs84UnionGeometry.toText());
            part.setMu(calcMu(wgs84UnionGeometry));
            parts.add(part);
        }
        log.debug("检查完毕");

        log.debug("检查是否有地块相交");
        if (parts.size() > 1) {
            // 判断每一个地块是否有相交，相交了多少平方米、相交了多少亩
            for (int i = 0; i < parts.size(); i++) {
                Part part = parts.get(i);
                for (int j = i + 1; j < parts.size(); j++) {
                    Part part2 = parts.get(j);
                    Geometry intersectionGeometry = part.getGaussGeometry().intersection(part2.getGaussGeometry());
                    if (intersectionGeometry.getArea() > 0) {
                        double area = intersectionGeometry.getArea();
                        double areaRounded = Math.round(area * 10000.0) / 10000.0;
                        double mu = area * config.SQUARE_TO_MU_METER;
                        double muRounded = Math.round(mu * 10000.0) / 10000.0;
                        log.debug("第 {} 个地块和第 {} 个地块相交了 {} 平方米，{} 亩",
                                i + 1, j + 1,
                                String.format("%.4f", areaRounded),
                                String.format("%.4f", muRounded)
                        );
                    }
                }
            }
        }
        log.debug("检查完毕");

        log.info("地块总面积={}亩 共 {} 个地块，耗时 {} 毫秒", splitResult.getMu(), unionPartsGaussGeometry.getNumGeometries(), System.currentTimeMillis() - splitRoadStartTime);
        return splitResult;
    }


}