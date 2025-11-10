package sunyu.util;

import java.time.Duration;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.concurrent.ConcurrentHashMap;
import java.util.stream.Collectors;

import org.geotools.geometry.jts.JTS;
import org.geotools.referencing.CRS;
import org.geotools.referencing.crs.DefaultGeographicCRS;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Envelope;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.LineString;
import org.locationtech.jts.geom.MultiPolygon;
import org.locationtech.jts.geom.Polygon;
import org.locationtech.jts.geom.prep.PreparedGeometry;
import org.locationtech.jts.geom.prep.PreparedGeometryFactory;
import org.opengis.referencing.crs.CoordinateReferenceSystem;
import org.opengis.referencing.operation.MathTransform;

import cn.hutool.log.Log;
import cn.hutool.log.LogFactory;
import sunyu.util.pojo.OutlinePart;
import sunyu.util.pojo.SplitRoadResult;
import sunyu.util.pojo.TrackPoint;

/**
 * GIS工具类，封装轨迹轮廓构建、道路分段、坐标转换、拓扑判断与面积计算。
 * 
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
        log.info("[构建{}] 开始", this.getClass().getSimpleName());
        log.info("[构建{}] 结束", this.getClass().getSimpleName());
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
        // 几何工厂
        private final GeometryFactory geometryFactory = new GeometryFactory();
        // 空几何集合 - 使用空数组参数更明确和安全
        private final Geometry EMPTYGEOM = geometryFactory.createGeometryCollection(new Geometry[0]);

        // WGS84坐标参考系统（全局常量，避免重复创建）
        private final CoordinateReferenceSystem WGS84_CRS = DefaultGeographicCRS.WGS84;

        // 高斯投影CRS缓存（线程安全），key格式："zone_falseEasting_centralMeridian"
        private final ConcurrentHashMap<String, CoordinateReferenceSystem> gaussCRSCache = new ConcurrentHashMap<>();

        // WGS84到高斯投影的坐标转换缓存（线程安全），key格式："zone_falseEasting_centralMeridian"
        private final ConcurrentHashMap<String, MathTransform> wgs84ToGaussTransformCache = new ConcurrentHashMap<>();

        // 高斯投影到WGS84的坐标转换缓存（线程安全），key格式："zone_falseEasting_centralMeridian"
        private final ConcurrentHashMap<String, MathTransform> gaussToWgs84TransformCache = new ConcurrentHashMap<>();

        // WGS84椭球长半轴（米），与Turf.js保持一致
        private final double EARTH_RADIUS = 6378137.0;

        // 最大速度（单位：米/秒），我们认为农机在田间作业，最大速度不会超过20米/秒
        private final double MAX_SPEED = 20.0;
    }

    /**
     * 构建器类，用于构建GisUtil实例。
     * <p>
     * 实现构建器模式，允许灵活配置和创建GisUtil对象，提供流式API进行参数设置。
     * </p>
     */
    public static class Builder {
        // 配置对象，包含各种常量和默认值
        private Config config = new Config();

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
     * 重置转换对象以释放资源。
     */
    @Override
    public void close() {
        log.info("[销毁{}] 开始", this.getClass().getSimpleName());
        log.info("[销毁{}] 结束", this.getClass().getSimpleName());
    }

    /**
     * 从缓存获取或创建高斯投影CRS
     * 
     * @param zone            投影带号
     * @param falseEasting    假东距
     * @param centralMeridian 中央经线
     * @return 高斯投影CRS对象，如果创建失败返回null
     */
    private CoordinateReferenceSystem getOrCreateGaussCRS(int zone, double falseEasting, double centralMeridian) {
        // 构建缓存key：zone_falseEasting_centralMeridian
        String cacheKey = String.format("%d_%.1f_%.1f", zone, falseEasting, centralMeridian);

        // 从缓存获取或创建高斯投影CRS
        return config.gaussCRSCache.computeIfAbsent(cacheKey, key -> {
            try {
                String wkt = String.format(
                        "PROJCS[\"WGS_1984_Gauss_Kruger_Zone_%d\", " +
                                "GEOGCS[\"GCS_WGS_1984\", " +
                                "DATUM[\"D_WGS_1984\", " +
                                "SPHEROID[\"WGS_1984\",6378137.0,298.257223563]], " +
                                "PRIMEM[\"Greenwich\",0.0], " +
                                "UNIT[\"Degree\",0.0174532925199433]], " +
                                "PROJECTION[\"Transverse_Mercator\"], " +
                                "PARAMETER[\"False_Easting\",%.1f], " +
                                "PARAMETER[\"False_Northing\",0.0], " +
                                "PARAMETER[\"Central_Meridian\",%.1f], " +
                                "PARAMETER[\"Scale_Factor\",1.0], " +
                                "PARAMETER[\"Latitude_Of_Origin\",0.0], " +
                                "UNIT[\"Meter\",1.0]]",
                        zone, falseEasting, centralMeridian);

                log.debug("高斯投影CRS WKT定义：{}", wkt);
                return CRS.parseWKT(wkt);
            } catch (Exception e) {
                log.warn("解析WKT失败：zone={}, falseEasting={}, centralMeridian={}, 错误={}",
                        zone, falseEasting, centralMeridian, e.getMessage());
                return null;
            }
        });
    }

    /**
     * 将WGS84坐标点转换为高斯投影坐标点
     * 支持全球范围，根据经度自动选择合适的高斯投影带
     * 
     * @param wgs84Point WGS84坐标系下的轨迹点
     * @return 高斯投影坐标系下的轨迹点
     */
    private TrackPoint wgs84PointTransformToGaussPoint(TrackPoint wgs84Point) {
        try {
            // 获取经度和纬度
            double longitude = wgs84Point.getLon();
            double latitude = wgs84Point.getLat();

            log.trace("=== 开始WGS84到高斯投影转换 ===");
            log.trace("原始坐标：经度={}, 纬度={}", longitude, latitude);

            // 计算高斯投影带号 (6度分带)
            // 全球范围: 经度-180到180，对应带号1-60
            int zone = (int) Math.floor((longitude + 180) / 6) + 1;
            log.trace("计算投影带号：zone={}", zone);

            // 计算中央经线
            double centralMeridian = (zone - 1) * 6 - 180 + 3;
            log.trace("计算中央经线：centralMeridian={}", centralMeridian);

            // 验证投影带号的合理性
            if (zone < 1 || zone > 60) {
                log.warn("投影带号超出合理范围：zone={}，经度={}", zone, longitude);
                return null;
            }

            // 创建高斯投影坐标参考系统（使用缓存避免重复解析WKT）
            // 全球支持模式：假东距包含带号信息，便于识别投影带
            // 格式：zone × 1000000 + 500000（如49带 = 49500000米）
            double falseEasting = zone * 1000000.0 + 500000.0;
            log.trace("计算假东距：falseEasting={}", falseEasting);

            // 从缓存获取或创建高斯投影CRS
            CoordinateReferenceSystem gaussCRS = getOrCreateGaussCRS(zone, falseEasting, centralMeridian);

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
            result.setSpeed(wgs84Point.getSpeed());
            result.setDirection(wgs84Point.getDirection());

            log.trace("=== WGS84到高斯投影转换完成 ===");
            return result;
        } catch (Exception e) {
            log.warn("坐标转换失败：经度={}, 纬度={}, 错误={}", wgs84Point.getLon(), wgs84Point.getLat(), e.getMessage());
            return null;
        }
    }

    private Geometry gaussGeometryToWgs84Geometry(Geometry gaussGeometry) {
        try {
            // 获取几何图形的边界信息来反推高斯投影参数
            Envelope env = gaussGeometry.getEnvelopeInternal();
            double centerX = (env.getMinX() + env.getMaxX()) / 2.0;

            log.debug("高斯投影几何到WGS84投影几何转换开始");
            log.trace("高斯几何边界：MinX={}, MaxX={}, MinY={}, MaxY={}",
                    env.getMinX(), env.getMaxX(), env.getMinY(), env.getMaxY());
            log.trace("几何中心X坐标：centerX={}", centerX);

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

            // 计算中央经线（与wgs84PointTransformToGaussPoint方法保持一致）
            double centralMeridian = (zone - 1) * 6 - 180 + 3;
            double falseEasting = zone * 1000000.0 + 500000.0;
            log.trace("计算投影参数：zone={}, centralMeridian={}, falseEasting={}", zone, centralMeridian, falseEasting);

            // 从缓存获取或创建高斯投影CRS
            CoordinateReferenceSystem gaussCRS = getOrCreateGaussCRS(zone, falseEasting, centralMeridian);

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
            log.debug("高斯投影几何到WGS84投影几何转换完成");
            return wgs84Geometry;
        } catch (Exception e) {
            log.warn("高斯投影几何到WGS84投影几何转换失败：错误={}", e.getMessage());
            return config.EMPTYGEOM;
        }
    }

    /**
     * 将高斯投影的几何图形转换为WGS84坐标系的WKT字符串
     * 
     * @param gaussGeometry 高斯投影坐标系下的几何图形
     * @return WGS84坐标系下的WKT字符串
     */
    private String gaussGeometryToWgs84WKT(Geometry gaussGeometry) {
        log.debug("高斯投影几何转换为WGS84投影WKT字符串开始");
        Geometry wgs84Geometry = gaussGeometryToWgs84Geometry(gaussGeometry);
        log.debug("高斯投影几何转换为WGS84投影WKT字符串完成");
        return wgs84Geometry.toText();
    }

    /**
     * 计算球面面积（平方米）
     * 使用球面多边形面积公式，与Turf.js的算法保持一致
     * 
     * @param geometry WGS84坐标系的几何图形
     * @return 球面面积（平方米）
     */
    private double calculateSphericalArea(Geometry geometry) {
        if (geometry == null || geometry.isEmpty()) {
            return 0.0;
        }

        double totalArea = 0.0;

        if (geometry instanceof Polygon) {
            totalArea = calculatePolygonSphericalArea((Polygon) geometry);
            log.debug("单多边形面积: {}平方米", totalArea);
        } else if (geometry instanceof MultiPolygon) {
            MultiPolygon multiPolygon = (MultiPolygon) geometry;
            for (int i = 0; i < multiPolygon.getNumGeometries(); i++) {
                Polygon polygon = (Polygon) multiPolygon.getGeometryN(i);
                double polyArea = calculatePolygonSphericalArea(polygon);
                totalArea += polyArea;
                log.debug("多边形{}面积: {}平方米", i, polyArea);
            }
            log.debug("MULTIPOLYGON总面积: {}平方米", totalArea);
        }

        return Math.abs(totalArea); // 确保面积为正值
    }

    /**
     * 计算单个多边形的球面面积（平方米）
     * 使用球面多边形面积公式，考虑地球曲率
     * 
     * @param polygon WGS84坐标系的多边形
     * @return 球面面积（平方米）
     */
    private double calculatePolygonSphericalArea(Polygon polygon) {
        // 外环面积
        double exteriorArea = calculateRingSphericalArea(polygon.getExteriorRing());
        log.trace("外环面积: {}平方米", exteriorArea);

        // 减去内环（孔洞）面积
        double holesArea = 0.0;
        for (int i = 0; i < polygon.getNumInteriorRing(); i++) {
            double holeArea = calculateRingSphericalArea(polygon.getInteriorRingN(i));
            holesArea += holeArea;
            log.trace("内环{}面积: {}平方米", i, holeArea);
        }

        double totalArea = exteriorArea - holesArea;
        log.trace("多边形总面积: {}平方米", totalArea);
        return totalArea;
    }

    /**
     * 计算线环的球面面积（平方米）
     * 使用球面多边形面积公式：A = R² * |Σ(λi+1 - λi) * sin(φi+1 + φi)/2|
     * 其中R为地球半径，λ为经度，φ为纬度
     * 
     * @param ring 线环（WGS84坐标系，单位：度）
     * @return 球面面积（平方米）
     */
    private double calculateRingSphericalArea(LineString ring) {
        if (ring == null || ring.isEmpty()) {
            return 0.0;
        }

        org.locationtech.jts.geom.Coordinate[] coords = ring.getCoordinates();
        if (coords.length < 3) {
            return 0.0; // 需要至少3个点才能形成多边形
        }

        double area = 0.0;

        // 使用正确的球面多边形面积公式（与Turf.js相同）
        for (int i = 0; i < coords.length - 1; i++) {
            double lon1 = Math.toRadians(coords[i].x);
            double lat1 = Math.toRadians(coords[i].y);
            double lon2 = Math.toRadians(coords[i + 1].x);
            double lat2 = Math.toRadians(coords[i + 1].y);

            // 正确的球面面积公式：A = R² × |Σ(λi+1 - λi) × sin((φi+1 + φi)/2)|
            area += (lon2 - lon1) * Math.sin((lat1 + lat2) / 2.0);
        }

        area = Math.abs(area) * config.EARTH_RADIUS * config.EARTH_RADIUS;
        return area;
    }

    /**
     * 直接从WKT字符串计算球面面积（平方米）
     * 使用与Turf.js相同的球面面积算法，支持POLYGON和MULTIPOLYGON
     * 
     * @param wkt WKT字符串（WGS84坐标系）
     * @return 球面面积（平方米）
     */
    private double calculateSphericalAreaFromWKT(String wkt) {
        if (wkt == null || wkt.trim().isEmpty()) {
            return 0.0;
        }

        wkt = wkt.trim();
        double totalArea = 0.0;

        try {
            if (wkt.startsWith("POLYGON")) {
                // 处理单个POLYGON
                totalArea = calculatePolygonAreaFromWKT(wkt);
                log.debug("POLYGON面积: {}平方米", totalArea);
            } else if (wkt.startsWith("MULTIPOLYGON")) {
                // 处理MULTIPOLYGON，提取所有多边形
                List<String> polygons = extractPolygonsFromMultiWKT(wkt);
                log.debug("MULTIPOLYGON包含 {} 个多边形", polygons.size());

                for (int i = 0; i < polygons.size(); i++) {
                    double polyArea = calculatePolygonAreaFromWKT(polygons.get(i));
                    totalArea += polyArea;
                    log.debug("多边形{}面积: {}平方米", i, polyArea);
                }
                log.debug("MULTIPOLYGON总面积: {}平方米", totalArea);
            }

        } catch (Exception e) {
            log.warn("WKT面积计算失败: {}", e.getMessage());
            return 0.0;
        }

        return Math.abs(totalArea);
    }

    /**
     * 从MULTIPOLYGON WKT中提取所有多边形WKT
     */
    private List<String> extractPolygonsFromMultiWKT(String multiWKT) {
        List<String> polygons = new ArrayList<>();

        // 移除MULTIPOLYGON前缀和后缀
        String content = multiWKT.substring(multiWKT.indexOf("(") + 1, multiWKT.lastIndexOf(")"));

        // 解析嵌套的多边形结构
        int level = 0;
        int start = 0;

        for (int i = 0; i < content.length(); i++) {
            char c = content.charAt(i);
            if (c == '(') {
                level++;
                if (level == 1) {
                    start = i;
                }
            } else if (c == ')') {
                level--;
                if (level == 0 && start >= 0) {
                    String polygonWKT = "POLYGON" + content.substring(start, i + 1);
                    polygons.add(polygonWKT);
                    start = -1;
                }
            }
        }

        return polygons;
    }

    /**
     * 计算单个POLYGON的面积（平方米）
     */
    private double calculatePolygonAreaFromWKT(String polygonWKT) {
        // 提取外环和内环坐标
        List<List<double[]>> rings = extractRingsFromPolygonWKT(polygonWKT);

        if (rings.isEmpty()) {
            return 0.0;
        }

        // 计算外环面积
        double exteriorArea = calculateRingArea(rings.get(0));

        // 减去内环（孔洞）面积
        double holesArea = 0.0;
        for (int i = 1; i < rings.size(); i++) {
            double holeArea = calculateRingArea(rings.get(i));
            holesArea += holeArea;
        }

        return exteriorArea - holesArea;
    }

    /**
     * 从POLYGON WKT中提取所有环（外环和内环）
     */
    private List<List<double[]>> extractRingsFromPolygonWKT(String polygonWKT) {
        List<List<double[]>> rings = new ArrayList<>();

        // 提取括号内的内容
        String content = polygonWKT.substring(polygonWKT.indexOf("(") + 1, polygonWKT.lastIndexOf(")"));

        // 解析环结构
        int level = 0;
        int start = 0;

        for (int i = 0; i < content.length(); i++) {
            char c = content.charAt(i);
            if (c == '(') {
                level++;
                if (level == 1) {
                    start = i + 1;
                }
            } else if (c == ')') {
                level--;
                if (level == 0 && start > 0) {
                    // 提取坐标字符串
                    String coordStr = content.substring(start, i);
                    List<double[]> coords = parseCoordinates(coordStr);
                    if (!coords.isEmpty()) {
                        rings.add(coords);
                    }
                    start = -1;
                }
            }
        }

        return rings;
    }

    /**
     * 解析坐标字符串为坐标列表
     */
    private List<double[]> parseCoordinates(String coordStr) {
        List<double[]> coordinates = new ArrayList<>();

        String[] coordPairs = coordStr.split(",");
        for (String pair : coordPairs) {
            pair = pair.trim();
            if (pair.isEmpty())
                continue;

            // 处理可能的多重空格
            String[] coords = pair.trim().split("\\s+");
            if (coords.length >= 2) {
                try {
                    double lon = Double.parseDouble(coords[0].trim());
                    double lat = Double.parseDouble(coords[1].trim());

                    // 验证坐标范围
                    if (lon >= -180 && lon <= 180 && lat >= -90 && lat <= 90) {
                        coordinates.add(new double[] { lon, lat });
                    } else {
                        log.warn("跳过无效坐标: 经度={}, 纬度={}", lon, lat);
                    }
                } catch (NumberFormatException e) {
                    log.warn("解析坐标失败: {}", pair);
                    continue;
                }
            }
        }

        // 确保多边形闭合
        if (coordinates.size() > 0) {
            double[] firstCoord = coordinates.get(0);
            double[] lastCoord = coordinates.get(coordinates.size() - 1);
            if (firstCoord[0] != lastCoord[0] || firstCoord[1] != lastCoord[1]) {
                coordinates.add(new double[] { firstCoord[0], firstCoord[1] });
            }
        }

        return coordinates;
    }

    /**
     * 计算单个环的面积（平方米）
     * 使用与Turf.js相同的球面面积算法
     */
    private double calculateRingArea(List<double[]> coordinates) {
        if (coordinates.size() < 3) {
            return 0.0;
        }

        double area = 0.0;

        // 使用正确的球面多边形面积公式：A = R² × |Σ(λi+1 - λi) × sin((φi+1 + φi)/2)|
        for (int i = 0; i < coordinates.size() - 1; i++) {
            double lon1 = Math.toRadians(coordinates.get(i)[0]);
            double lat1 = Math.toRadians(coordinates.get(i)[1]);
            double lon2 = Math.toRadians(coordinates.get(i + 1)[0]);
            double lat2 = Math.toRadians(coordinates.get(i + 1)[1]);

            area += (lon2 - lon1) * Math.sin((lat1 + lat2) / 2.0);
        }

        // 应用地球半径
        area = Math.abs(area) * config.EARTH_RADIUS * config.EARTH_RADIUS;

        log.debug("环面积: {}平方米, 坐标点数: {}", area, coordinates.size());
        return area;
    }

    public double calcMuByWgs84Geometry(Geometry wgs84Geometry) {
        if (wgs84Geometry == null || wgs84Geometry.isEmpty()) {
            return 0.0;
        }

        try {
            // 使用球面面积计算算法，与Turf.js对齐
            double areaSqm = calculateSphericalArea(wgs84Geometry);

            // 转换为亩：1亩 = 2000/3平方米，四舍五入保留4位小数
            return Math.round((areaSqm / (2000.0 / 3.0)) * 10000.0) / 10000.0;

        } catch (Exception e) {
            log.warn("WGS84几何图形计算亩数失败: {}", e.getMessage());
            return 0.0;
        }
    }

    public double calcMuByWgs84WKT(String wgs84WKT) {
        if (wgs84WKT == null || wgs84WKT.trim().isEmpty()) {
            log.warn("WKT字符串为空或null");
            return 0.0;
        }

        try {
            // 直接从WKT字符串计算球面面积，不转换为几何图形
            double areaSqm = calculateSphericalAreaFromWKT(wgs84WKT);

            // 转换为亩：1亩 = 2000/3平方米，四舍五入保留4位小数
            double mu = Math.round((areaSqm / (2000.0 / 3.0)) * 10000.0) / 10000.0;

            log.info("WKT面积计算结果: {}平方米 = {}亩", areaSqm, mu);
            return mu;

        } catch (Exception e) {
            log.warn("WKT字符串计算亩数失败: {}", e.getMessage());
            return 0.0;
        }
    }

    public SplitRoadResult splitRoad(List<TrackPoint> wgs84Points, double totalWidthM) {
        SplitRoadResult result = new SplitRoadResult();
        result.setTotalWidthM(totalWidthM);
        result.setOutline(config.EMPTYGEOM);
        result.setWkt(config.EMPTYGEOM.toText());

        log.debug("参数信息：轨迹点数量={} 总宽度={}米", wgs84Points.size(), totalWidthM);
        log.debug("过滤时间为空、经纬度异常、经纬度为0、速度为0、方向角异常的轨迹点");
        List<TrackPoint> filteredWgs84Points = new ArrayList<>();
        for (TrackPoint point : wgs84Points) {
            if (point.getTime() != null// 定位时间不能为空
                    && point.getLon() >= -180 && point.getLon() <= 180// 经度必须在[-180,180]之间
                    && point.getLat() >= -90 && point.getLat() <= 90// 纬度必须在[-90,90]之间
                    && point.getLon() != 0 && point.getLat() != 0// 经度和纬度不能为0
                    && point.getSpeed() > 0// 速度必须大于0
                    && point.getDirection() >= 0 && point.getDirection() <= 360) {// 方向角必须在[0,360]之间
                filteredWgs84Points.add(point);
            } else {
                log.trace("被过滤点位信息，时间：{} 经度：{} 纬度：{} 速度：{} 方向角：{}", point.getTime(), point.getLon(), point.getLat(),
                        point.getSpeed(), point.getDirection());
            }
        }
        log.debug("过滤后的轨迹点数量：{}", filteredWgs84Points.size());
        if (filteredWgs84Points.size() < 6) {
            log.warn("过滤后的轨迹点数量小于6个，无法进行拆分");
            return result;
        }

        log.debug("按定位时间升序排序");
        filteredWgs84Points.sort((p1, p2) -> p1.getTime().compareTo(p2.getTime()));

        // 获取上报时间间隔分布
        Map<Integer, Integer> intervalDistribution = new HashMap<>();
        log.info("统计轨迹点时间间隔分布，总点数：{}", filteredWgs84Points.size());
        for (int i = 1; i < filteredWgs84Points.size(); i++) {
            TrackPoint prevPoint = filteredWgs84Points.get(i - 1);
            TrackPoint currPoint = filteredWgs84Points.get(i);

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
        log.info("最小有效时间间隔：{}秒（点数最多的时间间隔）", minEffectiveInterval);

        int timeCutThreshold = minEffectiveInterval + 1; // 时间切割阈值
        log.debug("按照速度和时间间隔切分成多段轨迹，速度超过{}米/秒或者时间间隔超过{}秒，就拆分轨迹段", config.MAX_SPEED, timeCutThreshold);
        List<List<TrackPoint>> wgs84PointsSegments = new ArrayList<>();
        List<TrackPoint> currentSegment = new ArrayList<>();
        TrackPoint prevPoint = null;

        for (TrackPoint point : filteredWgs84Points) {
            boolean shouldSplit = false;

            // 速度超限切割
            if (point.getSpeed() >= config.MAX_SPEED) {
                shouldSplit = true;
            }

            // 时间间隔超限切割
            if (prevPoint != null) {
                Duration duration = Duration.between(prevPoint.getTime(), point.getTime());
                int timeDiffSeconds = (int) duration.getSeconds();
                if (timeDiffSeconds > timeCutThreshold) {
                    shouldSplit = true;
                }
            }

            if (shouldSplit) {
                if (!currentSegment.isEmpty()) {
                    wgs84PointsSegments.add(currentSegment);
                    currentSegment = new ArrayList<>();
                }
            }

            currentSegment.add(point);
            prevPoint = point;
        }
        if (!currentSegment.isEmpty()) {
            wgs84PointsSegments.add(currentSegment);
        }
        log.debug("按照速度和时间间隔切分后的轨迹段数量：{}", wgs84PointsSegments.size());

        log.debug("过滤掉段内小于6个点的轨迹段");
        wgs84PointsSegments.removeIf(segment -> segment.size() < 6);
        log.debug("过滤后的轨迹段数量：{}", wgs84PointsSegments.size());

        List<Geometry> gaussBufferedGeometries = new ArrayList<>();
        for (List<TrackPoint> wgs84PointsSegment : wgs84PointsSegments) {
            if (wgs84PointsSegment.size() < 6) {
                log.warn("轨迹段点数量小于6个，无法进行拆分");
                continue;
            }

            log.debug("将轨迹段中的WGS84坐标转换为高斯投影坐标 轨迹段点数量：{}", wgs84PointsSegment.size());
            List<TrackPoint> gaussPoints = new ArrayList<>();
            for (TrackPoint point : wgs84PointsSegment) {
                TrackPoint gaussPoint = wgs84PointTransformToGaussPoint(point);
                if (gaussPoint != null) {
                    gaussPoints.add(gaussPoint);
                }
            }
            log.debug("转换后的轨迹段点数量：{}", gaussPoints.size());
            if (gaussPoints.size() < 6) {
                log.warn("转换后的轨迹段点数量小于6个，无法进行拆分");
                continue;
            }

            log.debug("使用高斯投影的轨迹点进行线缓冲，左右缓冲总宽度：{}米", totalWidthM);
            try {
                // 1. 创建线几何：将高斯投影点转换为Coordinate数组
                Coordinate[] gaussCoordinates = new Coordinate[gaussPoints.size()];
                for (int i = 0; i < gaussPoints.size(); i++) {
                    TrackPoint point = gaussPoints.get(i);
                    gaussCoordinates[i] = new Coordinate(point.getLon(), point.getLat());
                }

                // 2. 创建LineString
                LineString lineString = config.geometryFactory.createLineString(gaussCoordinates);
                log.debug("创建线几何成功，点数：{}", gaussCoordinates.length);

                // 3. 执行缓冲操作：总宽度的一半作为缓冲距离
                double bufferDistance = totalWidthM / 2.0;
                Geometry gaussBufferedGeometry = lineString.buffer(bufferDistance);
                log.debug("线缓冲成功，缓冲距离：{}米，结果几何类型：{}", bufferDistance, gaussBufferedGeometry.getGeometryType());

                gaussBufferedGeometries.add(gaussBufferedGeometry);
            } catch (Exception e) {
                log.warn("线缓冲失败：总宽度={}米，错误={}", totalWidthM, e.getMessage());
            }
        }

        if (gaussBufferedGeometries.isEmpty()) {
            log.warn("没有成功创建任何区块，无法进行合并");
            return result;
        }

        log.debug("合并所有高斯投影缓冲几何");
        Geometry gaussUnionGeometry = config.geometryFactory
                .createGeometryCollection(gaussBufferedGeometries.toArray(new Geometry[0]))
                .union();
        log.debug("合并后的几何类型：{}，几何数量：{}", gaussUnionGeometry.getGeometryType(), gaussUnionGeometry.getNumGeometries());

        List<OutlinePart> outlineParts = new ArrayList<>();
        if (gaussUnionGeometry instanceof Polygon) {
            Geometry wgs84Geometry = gaussGeometryToWgs84Geometry(gaussUnionGeometry);

            // 筛选在合并几何图形内的WGS84点位 - 使用PreparedGeometry优化性能
            List<TrackPoint> wgs84PointsSegment = new ArrayList<>();

            // 获取几何图形的边界框，避免重复计算
            Envelope geometryEnvelope = wgs84Geometry.getEnvelopeInternal();

            // 使用PreparedGeometry预优化，大幅提升空间判断性能
            PreparedGeometry preparedWgs84Geometry = PreparedGeometryFactory.prepare(wgs84Geometry);

            // 记录性能日志
            long startTime = System.currentTimeMillis();

            // 使用并行流处理，结合PreparedGeometry提高大数据量场景下的性能
            wgs84PointsSegment = filteredWgs84Points.parallelStream()
                    .filter(point -> {
                        try {
                            // 第一步：边界框快速过滤 - 数值比较，性能极高
                            if (!geometryEnvelope.contains(point.getLon(), point.getLat())) {
                                return false;
                            }

                            // 第二步：使用PreparedGeometry进行精确空间判断 - 性能比直接contains提升10-100倍
                            return preparedWgs84Geometry.contains(config.geometryFactory.createPoint(
                                    new Coordinate(point.getLon(), point.getLat())));
                        } catch (Exception e) {
                            log.trace("点位空间判断失败：经度{} 纬度{} 错误：{}",
                                    point.getLon(), point.getLat(), e.getMessage());
                            return false;
                        }
                    })
                    .collect(Collectors.toList());

            long endTime = System.currentTimeMillis();
            log.debug("点位空间判断完成，原始点位数：{}，筛选后点位数：{}，耗时：{}ms",
                    filteredWgs84Points.size(), wgs84PointsSegment.size(), (endTime - startTime));

            OutlinePart outlinePart = new OutlinePart();
            outlinePart.setTotalWidthM(totalWidthM);
            outlinePart.setOutline(gaussUnionGeometry);
            outlinePart.setWkt(wgs84Geometry.toText());
            outlinePart.setMu(calcMuByWgs84Geometry(wgs84Geometry));
            outlinePart.setTrackPoints(wgs84PointsSegment);
            outlinePart.setStartTime(wgs84PointsSegment.get(0).getTime());
            outlinePart.setEndTime(wgs84PointsSegment.get(wgs84PointsSegment.size() - 1).getTime());
            outlineParts.add(outlinePart);
        } else if (gaussUnionGeometry instanceof MultiPolygon) {
            MultiPolygon multiPolygon = (MultiPolygon) gaussUnionGeometry;
            for (int i = 0; i < multiPolygon.getNumGeometries(); i++) {
                Polygon gaussGeometry = (Polygon) multiPolygon.getGeometryN(i);

                Geometry wgs84Geometry = gaussGeometryToWgs84Geometry(gaussGeometry);

                // 筛选在合并几何图形内的WGS84点位 - 使用PreparedGeometry优化性能
                List<TrackPoint> wgs84PointsSegment = new ArrayList<>();

                // 获取几何图形的边界框，避免重复计算
                Envelope geometryEnvelope = wgs84Geometry.getEnvelopeInternal();

                // 使用PreparedGeometry预优化，大幅提升空间判断性能
                PreparedGeometry preparedWgs84Geometry = PreparedGeometryFactory.prepare(wgs84Geometry);

                // 记录性能日志
                long startTime = System.currentTimeMillis();

                // 使用并行流处理，结合PreparedGeometry提高大数据量场景下的性能
                wgs84PointsSegment = filteredWgs84Points.parallelStream()
                        .filter(point -> {
                            try {
                                // 第一步：边界框快速过滤 - 数值比较，性能极大
                                if (!geometryEnvelope.contains(point.getLon(), point.getLat())) {
                                    return false;
                                }

                                // 第二步：使用PreparedGeometry进行精确空间判断 - 性能比直接contains提升10-100倍
                                return preparedWgs84Geometry.contains(config.geometryFactory.createPoint(
                                        new Coordinate(point.getLon(), point.getLat())));
                            } catch (Exception e) {
                                log.trace("点位空间判断失败：经度{} 纬度{} 错误：{}",
                                        point.getLon(), point.getLat(), e.getMessage());
                                return false;
                            }
                        })
                        .collect(Collectors.toList());

                long endTime = System.currentTimeMillis();
                log.debug("多边形{}点位空间判断完成，原始点位数：{}，筛选后点位数：{}，耗时：{}ms",
                        i, filteredWgs84Points.size(), wgs84PointsSegment.size(), (endTime - startTime));

                OutlinePart outlinePart = new OutlinePart();
                outlinePart.setTotalWidthM(totalWidthM);
                outlinePart.setOutline(gaussGeometry);
                outlinePart.setWkt(wgs84Geometry.toText());
                outlinePart.setMu(calcMuByWgs84Geometry(wgs84Geometry));
                outlinePart.setTrackPoints(wgs84PointsSegment);
                outlinePart.setStartTime(wgs84PointsSegment.get(0).getTime());
                outlinePart.setEndTime(wgs84PointsSegment.get(wgs84PointsSegment.size() - 1).getTime());
                outlineParts.add(outlinePart);
            }
        }

        // 4. 设置结果
        result.setOutline(gaussUnionGeometry);
        result.setWkt(gaussGeometryToWgs84WKT(gaussUnionGeometry));
        result.setParts(outlineParts);
        result.setMu(outlineParts != null ? outlineParts.stream()
                .filter(Objects::nonNull)
                .mapToDouble(OutlinePart::getMu)
                .sum() : 0.0);

        log.debug("最终生成区块数量：{}块，总亩数：{}亩", outlineParts.size(), result.getMu());

        return result;
    }

}