package sunyu.util;

import cn.hutool.core.collection.CollUtil;
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
import elki.distance.minkowski.SquaredEuclideanDistance;
import org.geotools.geometry.jts.JTS;
import org.geotools.geometry.jts.JTSFactoryFinder;
import org.geotools.referencing.CRS;
import org.geotools.referencing.crs.DefaultGeographicCRS;
import org.locationtech.jts.geom.*;
import org.locationtech.jts.index.strtree.STRtree;
import org.locationtech.jts.simplify.DouglasPeuckerSimplifier;
import org.opengis.referencing.crs.CoordinateReferenceSystem;
import org.opengis.referencing.operation.MathTransform;
import sunyu.util.pojo.GaussPoint;
import sunyu.util.pojo.Wgs84Point;

import java.time.Duration;
import java.util.*;
import java.util.concurrent.ConcurrentHashMap;
import java.util.stream.Collectors;

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
        private final Geometry EMPTYGEOM = GEOMETRY_FACTORY.createGeometryCollection();

        /**
         * WGS84坐标参考系统，使用默认地理CRS，避免重复创建
         * <p>这是全球标准的地理坐标参考系统，所有输入输出都基于此坐标系</p>
         */
        private final CoordinateReferenceSystem WGS84_CRS = DefaultGeographicCRS.WGS84;

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
         * 最小作业幅宽阈值（米），低于此值认为参数无效，用于农业轨迹处理
         * <p>设置此阈值是为了过滤掉不合理的作业幅宽参数，确保算法稳定性</p>
         */
        private final double MIN_WORKING_WIDTH_M = 1.0;

        /**
         * 最大返回多边形数量
         * <p>限制返回结果的数量，避免在处理复杂轨迹时产生过多碎片化的几何图形</p>
         */
        private final int MAX_GEOMETRY = 10;

        /**
         * 两点间最大作业距离(米)
         * <p>对应18km/h的作业速度（5米/秒），用于判断轨迹点之间的合理性</p>
         */
        private final double MAX_WORK_DISTANCE_M = 18 / 3.6;

        /**
         * 几何图形膨胀收缩距离（米）
         * <p>
         * 用于几何图形的膨胀-收缩优化操作。
         * 该距离影响几何图形的平滑效果，较大值平滑效果更强，较小值保持更多原始细节。
         * 先正向膨胀再反向收缩，用于去除小孔洞、平滑边界等几何优化。
         * </p>
         */
        private final double BUFFER_SMOOTHING_DISTANCE_M = 1;
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

    private CoordinateReferenceSystem getGaussCRS(int zone, double falseEasting, double centralMeridian) {
        // 创建缓存键，使用带号、假东距和中央经线作为唯一标识
        String cacheKey = String.format("%d_%.1f_%.1f", zone, falseEasting, centralMeridian);

        // 使用computeIfAbsent实现线程安全的缓存机制，只在缓存未命中时创建新CRS
        return config.GAUSS_CRS_CACHE.computeIfAbsent(cacheKey, key -> {
            try {
                log.debug("创建高斯投影CRS：投影带号={}, 假东距={}, 中央经线={}", zone, falseEasting,
                        centralMeridian);
                // 定义高斯-克吕格投影坐标系 - 使用预构建的WKT模板
                // WKT格式定义了完整的坐标系统参数，包括基准面、椭球体、投影方式和参数
                String wktTemplate = "PROJCS[\"Gauss_Kruger_ZONE_" + zone + "\", GEOGCS[\"GCS_WGS_1984\", DATUM[\"WGS_1984\", SPHEROID[\"WGS_84\", 6378137.0, 298.257223563]], PRIMEM[\"Greenwich\", 0.0], UNIT[\"Degree\", 0.0174532925199433]], PROJECTION[\"Transverse_Mercator\"], PARAMETER[\"False_Easting\", %.6f], PARAMETER[\"False_Northing\", 0.0], PARAMETER[\"Central_Meridian\", %.6f], PARAMETER[\"Scale_Factor\", 1.0], PARAMETER[\"Latitude_Of_Origin\", 0.0], UNIT[\"Meter\", 1.0]]";
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

    private Geometry processChunk(List<GaussPoint> points, int startIndex, int chunkSize,
                                  int totalPoints, double bufferWidth) {
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

        log.debug("处理块: {}-{}, 点数: {}, buffer几何类型: {}, 耗时: {}ms",
                startIndex, end, chunkLength, chunkBuffer.getGeometryType(),
                System.currentTimeMillis() - chunkStartTime);

        return chunkBuffer;
    }

    private Geometry mergeGeometriesRecursively(List<Geometry> geometries) {
        // 基本情况：列表为空或只有一个元素
        if (geometries.isEmpty()) {
            return config.EMPTYGEOM;
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
            return config.EMPTYGEOM;
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

        return config.EMPTYGEOM;
    }

    private Geometry processLargeSegmentInChunks(List<GaussPoint> points, double bufferWidth) {
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
            if (!geom.isEmpty()) {
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


    public Geometry toWgs84Geometry(Geometry gaussGeometry) {
        try {
            // 开始高斯投影几何到WGS84投影几何转换
            // 获取几何图形的边界信息来反推高斯投影参数
            Envelope env = gaussGeometry.getEnvelopeInternal();
            double centerX = (env.getMinX() + env.getMaxX()) / 2.0;

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
            // 计算中央经线（与wgs84PointTransformToGaussPoint方法保持一致）
            double centralMeridian = (zone - 1) * 6 - 180 + 3;
            double falseEasting = zone * 1000000.0 + 500000.0;

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
            MathTransform transform = config.GAUSS_TO_WGS84_TRANSFORM_CACHE.computeIfAbsent(cacheKey, key -> {
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
            Geometry wgs84Geometry = JTS.transform(gaussGeometry, transform);

            // 验证转换后的WGS84坐标合理性
            Envelope wgs84Env = wgs84Geometry.getEnvelopeInternal();
            if (wgs84Env.getMinX() < -180 || wgs84Env.getMaxX() > 180 ||
                    wgs84Env.getMinY() < -90 || wgs84Env.getMaxY() > 90) {
                log.warn("转换后的WGS84坐标超出合理范围：MinLon={}, MaxLon={}, MinLat={}, MaxLat={}",
                        wgs84Env.getMinX(), wgs84Env.getMaxX(), wgs84Env.getMinY(), wgs84Env.getMaxY());
                return config.EMPTYGEOM;
            }
            // 高斯投影几何到WGS84投影几何转换完成
            return wgs84Geometry;
        } catch (Exception e) {
            log.warn("高斯投影几何到WGS84投影几何转换失败：错误={}", e.getMessage());
            return config.EMPTYGEOM;
        }
    }

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

                // 获取或创建坐标转换对象
                String cacheKey = String.format("%d_%.1f_%.1f", zone, falseEasting, centralMeridian);
                CoordinateReferenceSystem gaussCRS = getGaussCRS(zone, falseEasting, centralMeridian);
                if (gaussCRS == null) {
                    continue;
                }

                MathTransform transform = config.WGS84_TO_GAUSS_TRANSFORM_CACHE.computeIfAbsent(cacheKey, key -> {
                    try {
                        return CRS.findMathTransform(config.WGS84_CRS, gaussCRS, true);
                    } catch (Exception e) {
                        log.warn("创建坐标转换失败：zone={}, falseEasting={}, centralMeridian={}, 错误={}",
                                zone, falseEasting, centralMeridian, e.getMessage());
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

    public void splitRoad(List<Wgs84Point> wgs84Wgs84Points, double workingWidth) {
        if (CollUtil.isEmpty(wgs84Wgs84Points)) {
            log.error("作业轨迹点列表不能为空");
            return;
        }
        if (workingWidth < 1) {
            log.error("作业幅宽必须大于 0 米");
            return;
        }
        log.debug("参数 wgs84TrackPoints 大小：{} workingWidth：{}", wgs84Wgs84Points.size(), workingWidth);

        double halfWorkingWidth = workingWidth / 2.0;

        log.debug("准备计算上报时间间隔分布");
        int minEffectiveInterval = 1; // 默认值
        Map<Integer, Integer> intervalDistribution = new HashMap<>();
        for (int i = 1; i < wgs84Wgs84Points.size(); i++) {
            Wgs84Point prevPoint = wgs84Wgs84Points.get(i - 1);
            Wgs84Point currPoint = wgs84Wgs84Points.get(i);

            // 计算时间间隔（秒）- LocalDateTime使用Duration
            Duration duration = Duration.between(prevPoint.getGpsTime(), currPoint.getGpsTime());
            int timeDiffSeconds = (int) duration.getSeconds();
            // 统计每个间隔的出现次数
            intervalDistribution.put(timeDiffSeconds, intervalDistribution.getOrDefault(timeDiffSeconds, 0) + 1);
        }
        // 获取最小有效时间间隔
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
        log.debug("最小有效上报时间间隔 {} 秒", minEffectiveInterval);

        log.debug("准备过滤异常点位信息");
        wgs84Wgs84Points = wgs84Wgs84Points.stream()
                .filter(p -> {
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
                })
                .collect(Collectors.toList());
        log.debug("过滤异常点位信息完成，剩余点位数量：{}", wgs84Wgs84Points.size());
        if (wgs84Wgs84Points.size() < 30) {
            log.error("作业轨迹点列表必须包含至少 30 个有效点位");
            return;
        }

        log.debug("准备对轨迹点按时间升序排序");
        wgs84Wgs84Points.sort(Comparator.comparing(Wgs84Point::getGpsTime));
        log.debug("轨迹点按时间升序排序完成");

        log.debug("准备去重前后两个经纬度点完全一致的轨迹点");
        List<Wgs84Point> filteredPoints = new ArrayList<>();
        Wgs84Point pre = null;
        for (Wgs84Point wgs84Point : wgs84Wgs84Points) {
            if (pre == null) {
                pre = wgs84Point;
                filteredPoints.add(pre);
                continue;
            }
            if (wgs84Point.getLongitude() == pre.getLongitude() && wgs84Point.getLatitude() == pre.getLatitude()) {
                log.warn("定位时间: {} 轨迹点经纬度与前一个点 {} 重复，抛弃", wgs84Point.getGpsTime(), pre.getGpsTime());
                continue;
            }
            pre = wgs84Point;
            filteredPoints.add(pre);
        }
        wgs84Wgs84Points = filteredPoints;
        log.debug("去重前后两个经纬度点完全一致的轨迹点完成，剩余点位数量：{}", wgs84Wgs84Points.size());
        if (wgs84Wgs84Points.size() < 30) {
            log.error("作业轨迹点列表必须包含至少 30 个有效点位");
            return;
        }

        // 转换为高斯投影坐标
        List<GaussPoint> gaussPoints = toGaussPointList(wgs84Wgs84Points);

        log.debug("准备创建StaticArrayDatabase");
        double eps = 5.0; // 1秒间隔保持6米（略大于最大速度5米/秒），10秒间隔保持35米
        int minPts = 4;

        log.debug("提取坐标数组");
        double[][] coords = new double[gaussPoints.size()][2];
        for (int i = 0; i < gaussPoints.size(); i++) {
            coords[i][0] = gaussPoints.get(i).getGaussX();
            coords[i][1] = gaussPoints.get(i).getGaussY();
        }

        log.debug("创建数据库");
        Database db = new StaticArrayDatabase(new ArrayAdapterDatabaseConnection(coords), null);
        db.initialize();

        log.debug("DBSCAN");
        DBSCAN<DoubleVector> dbscan = new DBSCAN<>(
                SquaredEuclideanDistance.STATIC, // 使用平方欧几里得距离
                eps,
                minPts
        );

        log.debug("获取Relation对象并运行聚类");
        Relation<DoubleVector> relation = db.getRelation(TypeUtil.DOUBLE_VECTOR_FIELD);
        Clustering<Model> result = dbscan.run(relation);

        log.debug("映射结果");
        DBIDRange ids = (DBIDRange) relation.getDBIDs();
        List<List<GaussPoint>> clusters = new ArrayList<>();
        for (Cluster<Model> cluster : result.getAllClusters()) {
            List<GaussPoint> clusterPoints = new ArrayList<>();
            for (DBIDIter iter = cluster.getIDs().iter(); iter.valid(); iter.advance()) {
                int offset = ids.getOffset(iter);
                clusterPoints.add(gaussPoints.get(offset));
            }
            clusters.add(clusterPoints);
        }
        log.debug("空间密集聚类完成，共发现 {} 个聚类簇", clusters.size());
        List<Geometry> unionGaussGeometries = new ArrayList<>();
        for (List<GaussPoint> cluster : clusters) {
            log.debug("聚类簇包含 {} 个点", cluster.size());
            if (cluster.size() > 30) {
                log.debug("聚类后，按时间升序排序");
                cluster.sort(Comparator.comparing(GaussPoint::getGpsTime));
                log.debug("创建线缓冲，缓冲半径：{} 米", halfWorkingWidth);
                Coordinate[] coordinates = cluster.stream()
                        .map(point -> new Coordinate(point.getGaussX(), point.getGaussY()))
                        .toArray(Coordinate[]::new);
                if (coordinates.length > 500) {
                    LineString lineString = config.GEOMETRY_FACTORY.createLineString(coordinates);
                    Geometry simplifiedGeometry = DouglasPeuckerSimplifier.simplify(lineString, 0.1);//0.1米容差
                    Coordinate[] simplifiedCoords = simplifiedGeometry.getCoordinates();
                    if (simplifiedCoords.length > 500) {
                        Geometry gaussGeometry = processLargeSegmentInChunks(cluster, halfWorkingWidth);
                        unionGaussGeometries.add(gaussGeometry);
                    } else {
                        Geometry gaussGeometry = simplifiedGeometry.buffer(halfWorkingWidth);
                        unionGaussGeometries.add(gaussGeometry);
                    }
                } else {
                    LineString lineString = config.GEOMETRY_FACTORY.createLineString(coordinates);
                    Geometry gaussGeometry = lineString.buffer(halfWorkingWidth);
                    unionGaussGeometries.add(gaussGeometry);
                }
            }
        }
        Geometry unionGaussGeometry = config.GEOMETRY_FACTORY
                .createGeometryCollection(unionGaussGeometries.toArray(new Geometry[0]))
                .union()
                .buffer(config.BUFFER_SMOOTHING_DISTANCE_M).buffer(-config.BUFFER_SMOOTHING_DISTANCE_M);
        Geometry wgs84UnionGeometry = toWgs84Geometry(unionGaussGeometry);
        log.debug("合并后的几何图形：{}", wgs84UnionGeometry.toText());
    }

}