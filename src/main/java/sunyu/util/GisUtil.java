package sunyu.util;

import java.time.Duration;
import java.time.LocalDateTime;
import java.time.ZoneOffset;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.concurrent.ConcurrentHashMap;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import org.geotools.geometry.DirectPosition2D;
import org.geotools.geometry.jts.JTS;
import org.geotools.referencing.CRS;
import org.geotools.referencing.crs.DefaultGeographicCRS;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.CoordinateSequence;
import org.locationtech.jts.geom.Envelope;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.LineString;
import org.locationtech.jts.geom.LinearRing;
import org.locationtech.jts.geom.MultiLineString;
import org.locationtech.jts.geom.MultiPolygon;
import org.locationtech.jts.geom.Point;
import org.locationtech.jts.geom.Polygon;
import org.locationtech.jts.geom.TopologyException;
import org.locationtech.jts.geom.impl.CoordinateArraySequence;
import org.locationtech.jts.geom.prep.PreparedGeometry;
import org.locationtech.jts.geom.prep.PreparedGeometryFactory;
import org.locationtech.jts.io.ParseException;
import org.locationtech.jts.io.WKTReader;
import org.locationtech.jts.operation.buffer.BufferOp;
import org.locationtech.jts.operation.buffer.BufferParameters;
import org.locationtech.jts.operation.overlay.OverlayOp;
import org.locationtech.jts.operation.overlay.snap.SnapIfNeededOverlayOp;
import org.locationtech.jts.operation.union.UnaryUnionOp;
import org.locationtech.jts.simplify.DouglasPeuckerSimplifier;
import org.locationtech.jts.simplify.TopologyPreservingSimplifier;
import org.opengis.referencing.crs.CoordinateReferenceSystem;
import org.opengis.referencing.operation.MathTransform;

import cn.hutool.core.collection.CollUtil;
import cn.hutool.core.util.StrUtil;
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
     * 将实际时间间隔映射到标准频率（正负2秒容错）
     * 
     * @param actualInterval      实际间隔秒数
     * @param standardFrequencies 标准频率数组
     * @return 对应的标准频率，如果无法匹配返回-1
     */
    private static int mapToStandardFrequency(int actualInterval, int[] standardFrequencies) {
        for (int freq : standardFrequencies) {
            if (Math.abs(actualInterval - freq) <= 2) {
                return freq;
            }
        }
        return -1;
    }

    /**
     * 计算两点间的球面距离（米）
     */
    private double haversine(double lat1, double lon1, double lat2, double lon2) {
        final double R = 6371000; // 地球半径（米）
        double latDistance = Math.toRadians(lat2 - lat1);
        double lonDistance = Math.toRadians(lon2 - lon1);
        double a = Math.sin(latDistance / 2) * Math.sin(latDistance / 2)
                + Math.cos(Math.toRadians(lat1)) * Math.cos(Math.toRadians(lat2))
                        * Math.sin(lonDistance / 2) * Math.sin(lonDistance / 2);
        double c = 2 * Math.atan2(Math.sqrt(a), Math.sqrt(1 - a));
        return R * c;
    }

    /**
     * 使用指定的breakDist构建线缓冲，支持自定义分段距离阈值
     */
    private Geometry buildSmartLineBufferWithBreakDist(List<TrackPoint> seg, double widthM, boolean enableLineBreak,
            double breakDist) throws Exception {
        if (seg == null || seg.size() < 3) {
            throw new IllegalArgumentException("轨迹段至少需要3个点");
        }

        // 坐标转换
        List<Coordinate> coords = new ArrayList<>();
        for (TrackPoint tp : seg) {
            coords.add(new Coordinate(tp.getLon(), tp.getLat()));
        }

        GeometryFactory gf = new GeometryFactory();
        CoordinateSequence coordSeq = new CoordinateArraySequence(coords.toArray(new Coordinate[0]));
        LineString line = new LineString(coordSeq, gf);

        // 使用指定的breakDist进行轨迹分段
        if (enableLineBreak && breakDist > 0) {
            double actualThreshold = breakDist * 2; // 提取阈值到变量，避免上下不一致
            log.debug("[buildSmartLineBufferWithBreakDist] 使用自定义breakDist={}m进行轨迹分段, 实际阈值={}m", breakDist,
                    actualThreshold);

            List<LineString> segments = new ArrayList<>();
            List<Coordinate> currentSegment = new ArrayList<>();
            currentSegment.add(coords.get(0));

            for (int i = 1; i < coords.size(); i++) {
                Coordinate prevCoord = coords.get(i - 1);
                Coordinate currCoord = coords.get(i);

                double distance = haversine(prevCoord.y, prevCoord.x, currCoord.y, currCoord.x);

                // 打印所有相邻点之间的距离用于调试
                log.trace("[buildSmartLineBufferWithBreakDist] 点{}到点{}距离={}m, 阈值={}m",
                        i - 1, i, distance, actualThreshold);

                // 如果当前点到下一点的距离超过阈值，则分段
                if (distance > actualThreshold) {
                    // 结束当前分段
                    if (currentSegment.size() >= 2) {
                        CoordinateSequence segmentSeq = new CoordinateArraySequence(
                                currentSegment.toArray(new Coordinate[0]));
                        segments.add(new LineString(segmentSeq, gf));
                        log.trace("[buildSmartLineBufferWithBreakDist] 创建分段: 点数={}, 距离超限={}m",
                                currentSegment.size(), distance);
                    }
                    // 开始新分段
                    currentSegment.clear();
                    currentSegment.add(prevCoord); // 包含前一个点
                    currentSegment.add(currCoord); // 包含当前点
                } else {
                    // 距离未超限，继续当前分段
                    if (currentSegment.isEmpty()) {
                        currentSegment.add(prevCoord);
                    }
                    currentSegment.add(currCoord);
                }
            }

            // 添加最后一个分段
            if (currentSegment.size() >= 2) {
                CoordinateSequence segmentSeq = new CoordinateArraySequence(currentSegment.toArray(new Coordinate[0]));
                segments.add(new LineString(segmentSeq, gf));
            }

            log.debug("[buildSmartLineBufferWithBreakDist] 轨迹分段完成: 分段数={}", segments.size());

            // 为每个分段构建缓冲
            log.debug("[buildSmartLineBufferWithBreakDist] 开始构建缓冲区: 分段数={}, 缓冲宽度={}m", segments.size(), widthM);
            long bufferStartTime = System.currentTimeMillis();
            List<Geometry> buffers = new ArrayList<>();
            for (int i = 0; i < segments.size(); i++) {
                LineString segment = segments.get(i);
                Geometry buffer = segment.buffer(widthM / 111000.0); // 转换为度
                buffers.add(buffer);
                log.trace("[buildSmartLineBufferWithBreakDist] 分段{}缓冲构建完成", i + 1);
            }
            long bufferEndTime = System.currentTimeMillis();
            log.debug("[buildSmartLineBufferWithBreakDist] 缓冲区构建完成: 耗时={}ms", bufferEndTime - bufferStartTime);

            // 合并所有缓冲（智能合并：只合并横向相交的缓冲区）
            log.debug("[buildSmartLineBufferWithBreakDist] 开始智能合并缓冲区: 缓冲区数量={}", buffers.size());
            long unionStartTime = System.currentTimeMillis();
            if (buffers.size() == 1) {
                long unionEndTime = System.currentTimeMillis();
                log.debug("[buildSmartLineBufferWithBreakDist] 缓冲区合并完成: 单个缓冲区无需合并, 耗时={}ms",
                        unionEndTime - unionStartTime);
                return buffers.get(0);
            } else {
                // 智能合并：只合并横向相交的缓冲区
                List<Geometry> mergedBuffers = smartMergeBuffers(buffers, segments);

                // 最终合并
                Geometry union = mergedBuffers.get(0);
                for (int i = 1; i < mergedBuffers.size(); i++) {
                    union = union.union(mergedBuffers.get(i));
                }
                long unionEndTime = System.currentTimeMillis();
                log.debug("[buildSmartLineBufferWithBreakDist] 缓冲区智能合并完成: 原始{}个缓冲区, 合并后{}个, 耗时={}ms",
                        buffers.size(), mergedBuffers.size(), unionEndTime - unionStartTime);
                return union;
            }
        } else {
            // 不分段，直接构建缓冲
            log.debug("[buildSmartLineBufferWithBreakDist] 不分段，直接构建缓冲");
            return line.buffer(widthM / 111000.0); // 转换为度
        }
    }

    /**
     * 智能合并缓冲区：只合并横向相交的缓冲区，不合并头尾相交的缓冲区
     * 
     * @param buffers  缓冲区列表
     * @param segments 对应的轨迹段列表
     * @return 合并后的缓冲区列表
     */
    private List<Geometry> smartMergeBuffers(List<Geometry> buffers, List<LineString> segments) {
        if (buffers.isEmpty()) {
            return new ArrayList<>();
        }

        List<Geometry> result = new ArrayList<>();
        boolean[] merged = new boolean[buffers.size()];

        for (int i = 0; i < buffers.size(); i++) {
            if (merged[i]) {
                continue;
            }

            Geometry currentBuffer = buffers.get(i);
            LineString currentSegment = segments.get(i);
            double currentDirection = calculateSegmentDirection(currentSegment);

            // 尝试与后续的缓冲区合并
            for (int j = i + 1; j < buffers.size(); j++) {
                if (merged[j]) {
                    continue;
                }

                // 检查是否相交
                if (!currentBuffer.intersects(buffers.get(j))) {
                    continue;
                }

                // 计算方向夹角
                LineString targetSegment = segments.get(j);
                double targetDirection = calculateSegmentDirection(targetSegment);
                double angleDiff = Math.abs(currentDirection - targetDirection);

                // 标准化角度差到0-90度范围
                if (angleDiff > 90) {
                    angleDiff = 180 - angleDiff;
                }

                // 判断是否为横向相交（夹角在60-120度之间）
                if (angleDiff >= 60 && angleDiff <= 120) {
                    // 横向相交，可以合并
                    try {
                        currentBuffer = currentBuffer.union(buffers.get(j));
                        merged[j] = true;
                        log.debug("[smartMergeBuffers] 合并横向相交的缓冲区 {} 和 {}，夹角={}度", i, j, angleDiff);
                    } catch (Exception e) {
                        log.warn("[smartMergeBuffers] 合并缓冲区 {} 和 {} 失败: {}", i, j, e.getMessage());
                    }
                } else {
                    log.debug("[smartMergeBuffers] 跳过头尾相交的缓冲区 {} 和 {}，夹角={}度", i, j, angleDiff);
                }
            }

            result.add(currentBuffer);
            merged[i] = true;
        }

        return result;
    }

    /**
     * 计算轨迹段的方向角度（相对于正北方向，单位：度）
     * 
     * @param segment 轨迹段
     * @return 方向角度（0-360度）
     */
    private double calculateSegmentDirection(LineString segment) {
        if (segment == null || segment.getNumPoints() < 2) {
            return 0.0;
        }

        Coordinate start = segment.getCoordinateN(0);
        Coordinate end = segment.getCoordinateN(segment.getNumPoints() - 1);

        // 计算相对于正北方向的角度
        double deltaX = end.x - start.x;
        double deltaY = end.y - start.y;

        double angle = Math.toDegrees(Math.atan2(deltaX, deltaY));

        // 标准化到0-360度
        if (angle < 0) {
            angle += 360;
        }

        return angle;
    }

    /**
     * 计算两个轨迹点之间的方向角度（相对于正北方向，单位：度）
     * 
     * @param start 起始轨迹点
     * @param end   结束轨迹点
     * @return 方向角度（0-360度）
     */
    private double calculateTrackPointDirection(TrackPoint start, TrackPoint end) {
        if (start == null || end == null) {
            log.warn("[calculateTrackPointDirection] 输入参数为空");
            return 0.0;
        }

        // 将WGS84坐标转换为米制坐标（简化版，使用第一个点作为基准）
        double deltaX = (end.getLon() - start.getLon()) * 111000.0 * Math.cos(Math.toRadians(start.getLat()));
        double deltaY = (end.getLat() - start.getLat()) * 111000.0;

        double angle = Math.toDegrees(Math.atan2(deltaX, deltaY));

        // 标准化到0-360度
        if (angle < 0) {
            angle += 360;
        }

        log.trace("[calculateTrackPointDirection] 起点({},{})->终点({},{})->角度={}度",
                start.getLon(), start.getLat(), end.getLon(), end.getLat(), angle);

        return angle;
    }

    /**
     * 跨窗口智能合并缓冲区：对所有窗口的最终结果再进行一次横向相交/相邻合并
     * 
     * 算法逻辑：
     * 1. 收集所有OutlinePart的几何图形和对应的轨迹点
     * 2. 使用与smartMergeBuffers相同的逻辑判断横向相交（60-120度夹角）
     * 3. 合并横向相交的缓冲区，保留独立的缓冲区
     * 
     * @param allParts 所有窗口的OutlinePart列表
     * @return 合并后的OutlinePart列表
     */
    private List<OutlinePart> smartMergeCrossWindowBuffers(List<OutlinePart> allParts) {
        if (allParts == null || allParts.isEmpty()) {
            log.debug("[smartMergeCrossWindowBuffers] 输入为空，直接返回");
            return new ArrayList<>();
        }
        if (allParts.size() == 1) {
            log.debug("[smartMergeCrossWindowBuffers] 只有一个区块，无需合并");
            return allParts;
        }

        log.debug("[smartMergeCrossWindowBuffers] 开始跨窗口合并，输入区块数={}", allParts.size());
        for (int i = 0; i < allParts.size(); i++) {
            OutlinePart part = allParts.get(i);
            log.debug("[smartMergeCrossWindowBuffers] 输入区块{}: 面积={}亩, 时间范围={}到{}",
                    i, part.getMu(), part.getStartTime(), part.getEndTime());
        }

        // 简单的一次性合并，避免死循环
        List<OutlinePart> result = new ArrayList<>();
        boolean[] merged = new boolean[allParts.size()];

        for (int i = 0; i < allParts.size(); i++) {
            if (merged[i])
                continue;

            OutlinePart currentPart = allParts.get(i);
            Geometry currentGeom = currentPart.getOutline();
            List<TrackPoint> currentPoints = currentPart.getTrackPoints();

            if (currentGeom == null || currentPoints == null || currentPoints.size() < 2) {
                result.add(currentPart);
                merged[i] = true;
                continue;
            }

            List<Geometry> buffersToMerge = new ArrayList<>();
            List<TrackPoint> allTrackPoints = new ArrayList<>(currentPoints);
            buffersToMerge.add(currentGeom);
            merged[i] = true;

            // 查找与当前区块相交的其他区块
            for (int j = i + 1; j < allParts.size(); j++) {
                if (merged[j])
                    continue;

                OutlinePart otherPart = allParts.get(j);
                Geometry otherGeom = otherPart.getOutline();
                List<TrackPoint> otherPoints = otherPart.getTrackPoints();

                if (otherGeom == null || otherPoints == null || otherPoints.size() < 2) {
                    continue;
                }

                // 只检查真正的几何相交或相邻（包括相切情况）
                boolean isIntersect = currentGeom.intersects(otherGeom) || currentGeom.touches(otherGeom);

                if (!isIntersect) {
                    log.debug("[smartMergeCrossWindowBuffers] 几何不相交，跳过合并");
                    continue;
                }

                log.debug("[smartMergeCrossWindowBuffers] 几何相交/相邻，当前区块面积={}亩, 其他区块面积={}亩",
                        calcMu(currentGeom), calcMu(otherGeom));
                buffersToMerge.add(otherGeom);
                allTrackPoints.addAll(otherPoints);
                merged[j] = true;
            }

            // 合并缓冲区
            if (buffersToMerge.size() > 1) {
                try {
                    log.debug("[smartMergeCrossWindowBuffers] 开始合并{}个缓冲区", buffersToMerge.size());
                    Geometry union = buffersToMerge.get(0);
                    for (int k = 1; k < buffersToMerge.size(); k++) {
                        union = union.union(buffersToMerge.get(k));
                        log.debug("[smartMergeCrossWindowBuffers] 合并第{}个缓冲区后，面积={}亩", k + 1, calcMu(union));
                    }

                    // 处理合并结果：根据几何类型创建相应的OutlinePart
                    allTrackPoints.sort(Comparator.comparing(TrackPoint::getTime));
                    LocalDateTime startTime = allTrackPoints.get(0).getTime();
                    LocalDateTime endTime = allTrackPoints.get(allTrackPoints.size() - 1).getTime();

                    if (union instanceof MultiPolygon) {
                        MultiPolygon multiPoly = (MultiPolygon) union;
                        log.debug("[smartMergeCrossWindowBuffers] 合并结果包含{}个多边形，为每个多边形创建OutlinePart",
                                multiPoly.getNumGeometries());

                        // 为MultiPolygon中的每个多边形创建单独的OutlinePart
                        for (int geomIndex = 0; geomIndex < multiPoly.getNumGeometries(); geomIndex++) {
                            Geometry geom = multiPoly.getGeometryN(geomIndex);
                            if (geom instanceof Polygon) {
                                Polygon poly = (Polygon) geom;
                                double polyMu = calcMu(poly);
                                String polyWkt = toWkt(poly);

                                // 构建并返回结果对象
                                OutlinePart polyPart = new OutlinePart();
                                polyPart.setOutline(poly);
                                polyPart.setStartTime(startTime);
                                polyPart.setEndTime(endTime);
                                polyPart.setMu(polyMu);
                                polyPart.setWkt(polyWkt);
                                polyPart.setTrackPoints(allTrackPoints);
                                polyPart.setTotalWidthM(currentPart.getTotalWidthM());
                                result.add(polyPart);
                                log.debug("[smartMergeCrossWindowBuffers] 创建多边形{}: 面积={}亩", geomIndex + 1, polyMu);
                            }
                        }
                    } else if (union instanceof Polygon) {
                        // 单个多边形的情况
                        double mergedMu = calcMu(union);
                        String mergedWkt = toWkt(union);

                        // 构建并返回结果对象
                        OutlinePart mergedPart = new OutlinePart();
                        mergedPart.setOutline(union);
                        mergedPart.setStartTime(startTime);
                        mergedPart.setEndTime(endTime);
                        mergedPart.setMu(mergedMu);
                        mergedPart.setWkt(mergedWkt);
                        mergedPart.setTrackPoints(allTrackPoints);
                        mergedPart.setTotalWidthM(currentPart.getTotalWidthM());
                        result.add(mergedPart);
                        log.debug("[smartMergeCrossWindowBuffers] 创建单个多边形: 面积={}亩", mergedMu);
                    } else {
                        log.warn("[smartMergeCrossWindowBuffers] 合并结果不是预期的多边形类型: {}", union.getGeometryType());
                    }
                    log.debug("[smartMergeCrossWindowBuffers] 合并了{}个横向相交的跨窗口区块", buffersToMerge.size());

                } catch (Exception e) {
                    log.warn("[smartMergeCrossWindowBuffers] 合并跨窗口缓冲区失败: {}", e.getMessage());
                    // 合并失败，添加原始区块
                    result.add(currentPart);
                }
            } else {
                // 没有需要合并的，直接添加原始区块
                log.debug("[smartMergeCrossWindowBuffers] 区块{}无需合并，直接添加，面积={}亩", i, calcMu(currentGeom));
                result.add(currentPart);
            }
        }

        log.debug("[smartMergeCrossWindowBuffers] 跨窗口合并完成，输出区块数={}", result.size());
        for (int i = 0; i < result.size(); i++) {
            OutlinePart part = result.get(i);
            log.debug("[smartMergeCrossWindowBuffers] 输出区块{}: 面积={}亩, 时间范围={}到{}",
                    i, part.getMu(), part.getStartTime(), part.getEndTime());
        }
        return result;
    }

    /**
     * 从几何体中提取所有多边形，支持MULTIPOLYGON和POLYGON类型
     */
    private List<Polygon> extractPolygons(Geometry geometry) {
        List<Polygon> polygons = new ArrayList<>();

        if (geometry instanceof Polygon) {
            polygons.add((Polygon) geometry);
        } else if (geometry instanceof MultiPolygon) {
            MultiPolygon multiPolygon = (MultiPolygon) geometry;
            for (int i = 0; i < multiPolygon.getNumGeometries(); i++) {
                Geometry subGeom = multiPolygon.getGeometryN(i);
                if (subGeom instanceof Polygon) {
                    polygons.add((Polygon) subGeom);
                }
            }
        }

        return polygons;
    }

    /**
     * 计算轨迹点列表中相邻点距离的中位数
     */
    private double calculateMedianAdjacentDistance(List<TrackPoint> points) {
        if (points.size() < 2) {
            return 0.0;
        }

        List<Double> distances = new ArrayList<>();
        for (int i = 1; i < points.size(); i++) {
            TrackPoint prev = points.get(i - 1);
            TrackPoint curr = points.get(i);
            double distance = haversine(prev.getLat(), prev.getLon(), curr.getLat(), curr.getLon());
            distances.add(distance);
        }

        Collections.sort(distances);

        int size = distances.size();
        if (size % 2 == 0) {
            return (distances.get(size / 2 - 1) + distances.get(size / 2)) / 2.0;
        } else {
            return distances.get(size / 2);
        }
    }

    /**
     * 智能窗口切分：基于时间和距离双重判断
     * 当相邻两点时间间隔>1分钟且距离>totalWidthM时切分窗口
     */
    private List<List<TrackPoint>> splitToTimeWindowsSmart(List<TrackPoint> points, double totalWidthM) {
        List<List<TrackPoint>> windows = new ArrayList<>();
        List<TrackPoint> currentWindow = new ArrayList<>();
        final long ONE_MINUTE_MS = 60 * 1000; // 1分钟

        log.debug("[splitToTimeWindowsSmart] 开始智能窗口切分: 总点数={}, 时间阈值={}分钟, 距离阈值={}m", points.size(),
                ONE_MINUTE_MS / 60000.0, totalWidthM);

        for (int i = 0; i < points.size(); i++) {
            TrackPoint point = points.get(i);

            if (currentWindow.isEmpty()) {
                currentWindow.add(point);
            } else {
                TrackPoint lastPoint = currentWindow.get(currentWindow.size() - 1);
                long timeDiff = point.getTime().toInstant(ZoneOffset.UTC).toEpochMilli()
                        - lastPoint.getTime().toInstant(ZoneOffset.UTC).toEpochMilli();
                double distance = haversine(lastPoint.getLat(), lastPoint.getLon(), point.getLat(), point.getLon());

                // 双重判断：时间>1分钟且距离>totalWidthM
                if (timeDiff > ONE_MINUTE_MS && distance > totalWidthM) {
                    log.debug("[splitToTimeWindowsSmart] 切分窗口: 时间间隔={}分钟, 距离={}m",
                            timeDiff / 60000.0, distance);

                    // 保存当前窗口
                    if (currentWindow.size() >= 3) {
                        windows.add(new ArrayList<>(currentWindow));
                    }
                    currentWindow.clear();
                }
                currentWindow.add(point);
            }
        }

        // 添加最后一个窗口
        if (currentWindow.size() >= 3) {
            windows.add(currentWindow);
        }

        return windows;
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

        // 预定义的空几何体，避免重复创建
        private final Geometry EMPTY_GEOMETRY = geometryFactory.createGeometryCollection(null);

        // 空几何体的WKT表示，用于初始化空结果
        private final String EMPTY_WKT = "GEOMETRYCOLLECTION EMPTY";

        // 默认轮廓返回的最多多边形数量（TopN）
        private final int DEFAULT_MAX_OUTLINE_SEGMENTS = 10;

        // 圆近似细分：让小面积轨迹缓冲区更圆滑（半圆效果）
        private final int DEFAULT_BUFFER_QUADRANT = 2;

        // 作业最大速度阈值（km/h），用于前置速度过滤；可通过 Builder 配置
        private final double WORK_MAX_SPEED_KMH = 18.0;
        // 作业最小速度阈值（km/h），用于前置速度过滤；可通过 Builder 配置
        private final double MIN_WORK_SPEED_KMH = 0.1;
        // 最小亩数动态阈值（亩），用于 splitRoad 动态过滤小块
        private final double MIN_MU_DYNAMIC_THRESHOLD_MU = 0.32;

        // 外缘细长裁剪开关（仅 splitRoad 使用）
        private final boolean ENABLE_OUTER_THIN_TRIM = true;
        // 拓扑简化容差系数（外缘细长裁剪用）：控制拓扑保留简化的精度，值越小精度越高
        private final double TOPOLOGY_SIMPLIFY_TOLERANCE_FACTOR = 0.2;
        // 拓扑简化容差上限（米）：防止大半径时容差过大
        private final double TOPOLOGY_SIMPLIFY_MAX_TOLERANCE_M = 1.0;
        // 外缘细长裁剪半径系数（相对单侧宽度）
        private final double THIN_TRIM_RADIUS_FACTOR = 1.5;
        // 宽幅到半径系数的映射表（使用TreeMap支持范围查找）
        private final TreeMap<Double, Double> WIDTH_TO_RADIUS_FACTOR_MAP = new TreeMap<Double, Double>() {
            {
                // 初始化默认映射关系
                put(1.0, 0.45);
                put(1.75, 1.4);
                put(2.5, 1.6);
                put(2.6, 1.8);
                put(2.7, 2.2);
                put(2.8, 2.4);
                put(3.0, 2.6);
            }
        };

        // 线段断裂控制开关（buildSmartLineBuffer 使用）
        private final boolean ENABLE_LINE_BREAK = true;
        // 线段断裂距离系数（倍数*单侧宽度），超过则切分会话
        private final double LINE_BREAK_FACTOR = 4;

        // 线简化控制开关（buildSmartLineBuffer 使用）
        private final boolean ENABLE_LINE_SIMPLIFY = true;
        // 线简化公差系数（倍数*单侧宽度），用于Douglas-Peucker
        private final double LINE_SIMPLIFY_TOL_FACTOR = 2;

        // 线缓冲样式（更圆滑边界）：拐角样式、端头样式与mitre限制
        private final int BUFFER_JOIN_STYLE = BufferParameters.JOIN_ROUND;
        // 线缓冲的端头样式（圆端），影响线段端点的缓冲形状
        private final int BUFFER_END_CAP_STYLE = BufferParameters.CAP_ROUND;
        // 线缓冲的mitre限制值，控制锐角的处理方式，值越小角越圆滑，值越大允许越尖的角
        private final double BUFFER_MITRE_LIMIT = 2.0;

        // 轨迹分段最小切割距离（米），用于splitRoad方法中的轨迹点分组
        private final TreeMap<Double, Double> MIN_SEGMENT_DISTANCE_THRESHOLD_MAP = new TreeMap<Double, Double>() {
            {
                put(1.0, 3.3);
                put(1.75, 2.7);
                put(2.8, 2.8);
                put(Double.MAX_VALUE, 4.0);
            }
        };

        // 多边形合并缓冲距离（米），用于调整空间合并的敏感度
        // 值越大，相邻多边形越容易合并；值越小，合并条件越严格
        private final double POLYGON_MERGE_BUFFER_DISTANCE_M = 0.5;

        // 兜底方案控制：当splitRoad结果为空或面积过小时启用getOutline重算
        private final boolean ENABLE_FALLBACK_TO_OUTLINE = true;
        // 兜底面积阈值（亩）：小于该值时触发兜底方案，使用getOutline重算
        private final double FALLBACK_MU_THRESHOLD = 1.5;
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
     * @param wgs84Lon WGS84经度（度），用于选择分带
     * @return 逆向投影变换（Gauss→WGS84）
     */
    private MathTransform getTxGaussToWgsByLon(double wgs84Lon) throws Exception {
        int zone = gkZoneFromLon(wgs84Lon);
        String key = "G2W:" + zone;
        MathTransform cached = config.txCache.get(key);
        if (cached != null)
            return cached;
        CoordinateReferenceSystem src = getGaussKrugerCRSByLon(wgs84Lon);
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

            // 构建并返回结果对象
            OutlinePart op = new OutlinePart();
            op.setOutline(polygon);
            op.setStartTime(startTime);
            op.setEndTime(endTime);
            op.setMu(mu);
            op.setWkt(wkt);
            op.setTrackPoints(pointsInPolygon);
            op.setTotalWidthM(totalWidthM);
            return op;
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
            SplitRoadResult result = new SplitRoadResult();
            result.setOutline(outline);
            result.setParts(parts);
            result.setWkt(outlineWkt);
            result.setTotalWidthM(totalWidthM);

            log.info("[splitRoad] 兜底方案执行成功: 区块数={}, 总亩数={}",
                    result.getParts().size(), result.getMu());

            return result;

        } catch (Exception e) {
            log.trace("[splitRoad] 兜底方案执行失败", e);

            // 兜底方案失败时，返回包含空几何的SplitRoadResult，避免调用者报错
            Geometry emptyGeometry = config.geometryFactory.createGeometryCollection(null);
            String emptyWkt = config.EMPTY_WKT;
            List<OutlinePart> emptyParts = new ArrayList<>();

            SplitRoadResult emptyResult = new SplitRoadResult();
            emptyResult.setOutline(emptyGeometry);
            emptyResult.setParts(emptyParts);
            emptyResult.setWkt(emptyWkt);
            emptyResult.setTotalWidthM(totalWidthM);

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
                return config.EMPTY_WKT;
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
            return g != null ? g.toText() : config.EMPTY_WKT;
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

        // 坐标/几何修复（自相交等）
        if (!hasValidCoordinates(g1)) {
            g1 = g1.buffer(0);
            if (g1.isEmpty()) {
                return new WktIntersectionResult(null, 0.0);
            }
        }
        if (!hasValidCoordinates(g2)) {
            g2 = g2.buffer(0);
            if (g2.isEmpty()) {
                return new WktIntersectionResult(null, 0.0);
            }
        }

        // 几何简化以减少精度问题
        double tolerance = 1e-8; // WGS84坐标系下的简化容差
        g1 = DouglasPeuckerSimplifier.simplify(g1, tolerance);
        g2 = DouglasPeuckerSimplifier.simplify(g2, tolerance);

        Geometry inter;
        try {
            // 使用SnapIfNeededOverlayOp处理拓扑异常
            SnapIfNeededOverlayOp snapOp = new SnapIfNeededOverlayOp(g1, g2);
            inter = snapOp.getResultGeometry(OverlayOp.INTERSECTION);
        } catch (TopologyException e) {
            log.warn("SnapIfNeededOverlayOp失败，尝试使用缓冲交集: {}", e.getMessage());
            // 如果SnapIfNeededOverlayOp失败，使用缓冲交集作为备选方案
            double bufferDistance = 1e-10; // 极小的缓冲距离
            Geometry g1Buffer = g1.buffer(bufferDistance);
            Geometry g2Buffer = g2.buffer(bufferDistance);
            inter = g1Buffer.intersection(g2Buffer);
        }

        if (inter == null || inter.isEmpty()) {
            return new WktIntersectionResult(null, 0.0);
        }
        String wkt = toWkt(inter);
        double mu = calcMuByWgs84Wkt(wkt);
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
        OutlinePart op = new OutlinePart();
        op.setOutline(outline);
        op.setStartTime(startTime);
        op.setEndTime(endTime);
        op.setMu(mu);
        op.setWkt(wkt);
        op.setTrackPoints(points);
        op.setTotalWidthM(totalWidthM);
        return op;
    }

    /**
     * 内部方法：执行splitRoad的核心逻辑
     * 
     * @param seg             轨迹点列表
     * @param totalWidthM     总宽度
     * @param maxSegments     最大区块数
     * @param enableLineBreak 是否启用线分割（多线程安全，使用临时参数）
     * @return 拆分结果
     * @throws Exception 处理异常
     */
    private SplitRoadResult doSplitRoad(List<TrackPoint> seg, double totalWidthM, Integer maxSegments,
            boolean enableLineBreak) throws Exception {
        long startTime = System.currentTimeMillis();
        double halfWidth = totalWidthM / 2.0;
        int segmentLimit = (maxSegments == null || maxSegments <= 0) ? config.DEFAULT_MAX_OUTLINE_SEGMENTS
                : maxSegments;

        log.debug("[splitRoad] 开始处理: 点数={}, 宽度={}m, 上限={}, enableLineBreak={}", seg.size(), totalWidthM, segmentLimit,
                enableLineBreak);

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
        Geometry outline = buildSmartLineBuffer(validPoints, halfWidth, enableLineBreak);
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

        SplitRoadResult result = new SplitRoadResult();
        result.setOutline(trimmed);
        result.setParts(parts);
        result.setWkt(outlineWkt);
        result.setTotalWidthM(totalWidthM);

        return result;
    }

    public SplitRoadResult splitRoad(List<TrackPoint> wgs84Points, double totalWidthM) {
        return splitRoad(wgs84Points, totalWidthM, config.DEFAULT_MAX_OUTLINE_SEGMENTS);
    }

    public SplitRoadResult splitRoad(List<TrackPoint> wgs84Points, double totalWidthM, int maxSegments) {
        SplitRoadResult result = new SplitRoadResult();
        result.setTotalWidthM(totalWidthM);
        result.setOutline(config.EMPTY_GEOMETRY);
        result.setWkt(config.EMPTY_WKT);
        result.setMu(0);

        // 参数校验，seg必须大于3个点
        if (CollUtil.isEmpty(wgs84Points) || wgs84Points.size() < 3) {
            log.warn("轨迹点列表必须包含至少3个有效点");
            return result;
        }
        // 参数校验，totalWidthM必须大于1米
        if (totalWidthM < 1) {
            log.warn("作业总宽度不能小于1米");
            return result;
        }
        // 参数校验，maxSegments必须大于0
        if (maxSegments <= 0) {
            log.warn("返回的区块数量上限必须大于0");
            return result;
        }

        log.debug("参数 wgs84Points.size={}, totalWidthM={}, maxSegments={}", wgs84Points.size(), totalWidthM,
                maxSegments);
        // 循环seg，如果时间为null，抛弃该点
        wgs84Points.removeIf(p -> p.getTime() == null);
        // 先将seg进行按时间升序排序
        wgs84Points.sort(Comparator.comparing(TrackPoint::getTime));
        // 过滤经纬度异常的点（纬度[-90,90]，经度[-180,180]）
        wgs84Points.removeIf(p -> p.getLat() < -90 || p.getLat() > 90 || p.getLon() < -180 || p.getLon() > 180);
        // 过滤经纬度为0的点
        wgs84Points.removeIf(p -> p.getLat() == 0 && p.getLon() == 0);
        if (CollUtil.isEmpty(wgs84Points) || wgs84Points.size() < 3) {
            log.warn("过滤后轨迹点列表不足3个点");
            return result;
        }
        log.debug("过滤后 wgs84Points.size={}", wgs84Points.size());

        // 将wgs84Points进行wgs84到高斯投影的转换
        List<TrackPoint> gaussPoints = convertWgs84ToGaussProjection(wgs84Points);
        if (CollUtil.isEmpty(gaussPoints) || gaussPoints.size() < 3) {
            log.warn("转换后的轨迹点列表不足3个有效点");
            return result;
        }

        // 过滤速度范围的点（高斯投影坐标系，米单位）
        /*
         * gaussPoints = filterBySpeedRangeMeter(gaussPoints, config.MIN_WORK_SPEED_KMH,
         * config.WORK_MAX_SPEED_KMH);
         * if (CollUtil.isEmpty(gaussPoints) || gaussPoints.size() < 3) {
         * log.warn("速度过滤后轨迹点列表不足3个有效点");
         * return result;
         * }
         */

        // 计算单侧缓冲宽度（totalWidthM是总宽度，需要除以2）
        double halfWidth = totalWidthM / 2.0;
        double halfHalfWidth = halfWidth / 2.0;

        // 实现按上报频率分组：统计速度过滤后的轨迹点能分出多少频率相同的组
        // gaussPoints已经在前面按时间排序（第2658行），无需再次排序
        if (gaussPoints.size() < 2) {
            log.warn("速度过滤后轨迹点数不足2个，无法进行频率分组");
            return result;
        }

        log.info("速度过滤后轨迹点数: {}", gaussPoints.size());

        // 第一步：定义标准频率模式（1秒、5秒、10秒、15秒）
        int[] standardFrequencies = { 1, 5, 10, 15 };
        Map<Integer, Integer> timeIntervalStats = new LinkedHashMap<>(); // 标准频率 -> 出现次数

        for (int i = 1; i < gaussPoints.size(); i++) {
            long timeDiff = java.time.Duration.between(gaussPoints.get(i - 1).getTime(), gaussPoints.get(i).getTime())
                    .getSeconds();

            // 只处理1-17秒的间隔（15秒+2秒容错）
            if (timeDiff >= 1 && timeDiff <= 17) {
                int standardFreq = mapToStandardFrequency((int) timeDiff, standardFrequencies);
                if (standardFreq != -1) {
                    timeIntervalStats.merge(standardFreq, 1, Integer::sum);
                }
            }
        }

        if (timeIntervalStats.isEmpty()) {
            log.warn("未找到符合标准频率的时间间隔");
            return result;
        }

        // 按标准频率排序输出
        List<Integer> sortedFrequencies = new ArrayList<>(timeIntervalStats.keySet());
        Collections.sort(sortedFrequencies);

        log.info("发现的标准频率模式（共{}种）:", sortedFrequencies.size());
        int totalIntervals = timeIntervalStats.values().stream().mapToInt(Integer::intValue).sum();

        for (int freq : sortedFrequencies) {
            int count = timeIntervalStats.get(freq);
            double percentage = (count * 100.0) / totalIntervals;
            log.info("频率{}秒: 出现{}次 ({}%)", freq, count, percentage);
        }

        // 第二步：根据发现的频率模式进行分组
        Map<Integer, List<List<TrackPoint>>> frequencyGroups = new LinkedHashMap<>();
        List<TrackPoint> currentGroup = new ArrayList<>();
        int currentFrequency = -1;

        for (int i = 0; i < gaussPoints.size(); i++) {
            TrackPoint currentPoint = gaussPoints.get(i);

            if (i == 0) {
                currentGroup.add(currentPoint);
                continue;
            }

            // 计算与前一个点的时间间隔（秒）
            long timeDiff = Duration.between(gaussPoints.get(i - 1).getTime(), currentPoint.getTime())
                    .getSeconds();

            // 如果时间间隔超过15秒，作为异常间隔处理（结束当前组）
            if (timeDiff > 15) {
                log.warn("发现时间间隔{}秒超过15秒，作为异常间隔处理", timeDiff);
                // 异常间隔，结束当前组
                if (currentGroup.size() >= 2) {
                    frequencyGroups.computeIfAbsent(currentFrequency, k -> new ArrayList<>())
                            .add(new ArrayList<>(currentGroup));
                }
                currentGroup.clear();
                currentGroup.add(currentPoint);
                currentFrequency = -1;
                continue;
            }

            // 将实际间隔映射到标准频率（正负2秒容错）
            int standardFreq = mapToStandardFrequency((int) timeDiff, standardFrequencies);

            if (standardFreq != -1) {
                if (currentFrequency == -1) {
                    currentFrequency = standardFreq;
                }

                if (standardFreq == currentFrequency) {
                    currentGroup.add(currentPoint);
                } else {
                    // 频率变化，结束当前组
                    if (currentGroup.size() >= 2) {
                        frequencyGroups.computeIfAbsent(currentFrequency, k -> new ArrayList<>())
                                .add(new ArrayList<>(currentGroup));
                    }
                    currentGroup.clear();
                    currentGroup.add(currentPoint);
                    currentFrequency = standardFreq;
                }
            } else {
                // 未识别的间隔（超过15秒或不在标准频率范围内），结束当前组并丢弃这组数据
                if (currentGroup.size() >= 2) {
                    frequencyGroups.computeIfAbsent(currentFrequency, k -> new ArrayList<>())
                            .add(new ArrayList<>(currentGroup));
                }
                currentGroup.clear();
                currentGroup.add(currentPoint);
                currentFrequency = -1;
            }
        }

        // 处理最后一组
        if (currentGroup.size() >= 2 && currentFrequency != -1) {
            frequencyGroups.computeIfAbsent(currentFrequency, k -> new ArrayList<>())
                    .add(new ArrayList<>(currentGroup));
        }

        int totalGroups = 0;
        int totalGroupPoints = 0;

        // 按标准频率顺序输出
        for (int freq : standardFrequencies) {
            if (!frequencyGroups.containsKey(freq)) {
                continue; // 跳过没有数据的频率
            }

            List<List<TrackPoint>> groups = frequencyGroups.get(freq);
            int groupCount = groups.size();
            int pointCount = groups.stream().mapToInt(List::size).sum();

            totalGroups += groupCount;
            totalGroupPoints += pointCount;

            log.info("频率{}秒: {}个组, {}个点", freq, groupCount, pointCount);

            // 打印每个组的详细信息
            for (int i = 0; i < groups.size(); i++) {
                List<TrackPoint> group = groups.get(i);
                if (!group.isEmpty()) {
                    log.debug("  组{}: {}个点, 时间:{} ~ {}",
                            i + 1,
                            group.size(),
                            group.get(0).getTime(),
                            group.get(group.size() - 1).getTime());
                }
            }
        }

        log.info("总频率模式数: {}", frequencyGroups.size());
        log.info("总组数: {}", totalGroups);

        // 获取最小间隔的分组
        List<Geometry> allPolygons = new ArrayList<>();
        if (!frequencyGroups.isEmpty()) {
            int minFrequency = frequencyGroups.keySet().stream().min(Integer::compareTo).orElse(0);
            List<List<TrackPoint>> minFreqGroups = frequencyGroups.get(minFrequency);
            int minFreqGroupCount = minFreqGroups.size();
            int minFreqPointCount = minFreqGroups.stream().mapToInt(List::size).sum();
            log.info("最小间隔{}秒: {}个组, {}个点", minFrequency, minFreqGroupCount, minFreqPointCount);

            for (int k = 0; k < minFreqGroups.size(); k++) {
                List<TrackPoint> group = minFreqGroups.get(k);
                // 计算每一组的距离中位数
                if (group.size() >= 3) {
                    List<Double> distances = new ArrayList<>();
                    for (int j = 1; j < group.size(); j++) {
                        TrackPoint prev = group.get(j - 1);
                        TrackPoint curr = group.get(j);
                        double distance = Math.hypot(
                                curr.getLon() - prev.getLon(),
                                curr.getLat() - prev.getLat());
                        distances.add(distance);
                    }

                    if (!distances.isEmpty()) {
                        distances.sort(Double::compareTo);
                        double medianDistance;
                        if (distances.size() % 2 == 0) {
                            medianDistance = (distances.get(distances.size() / 2 - 1) +
                                    distances.get(distances.size() / 2)) / 2.0;
                        } else {
                            medianDistance = distances.get(distances.size() / 2);
                        }

                        // 计算平均数、最大距离、最小距离
                        double sumDistance = distances.stream().mapToDouble(Double::doubleValue).sum();
                        double avgDistance = sumDistance / distances.size();
                        double maxDistance = distances.get(distances.size() - 1);

                        log.info("  组{}: {}个点, 距离统计: 中位数={}米, 平均数={}米, 最大={}米",
                                k + 1, group.size(), medianDistance, avgDistance, maxDistance);

                        if (0 < medianDistance && medianDistance <= 20) {
                            try {
                                double segmentDistanceThreshold = Math.max(
                                        config.MIN_SEGMENT_DISTANCE_THRESHOLD_MAP.floorEntry(totalWidthM).getValue(),
                                        minFrequency * totalWidthM);

                                List<List<TrackPoint>> segments = new ArrayList<>();
                                List<TrackPoint> currentSegment = new ArrayList<>();

                                log.debug("开始按{}米距离阈值对组{}进行轨迹分段", segmentDistanceThreshold, k + 1);

                                // 使用当前组的点进行分段，而不是所有轨迹点
                                for (int i = 0; i < group.size(); i++) {
                                    TrackPoint currentPoint = group.get(i);

                                    if (currentSegment.isEmpty()) {
                                        // 第一段或新段开始的第一个点
                                        currentSegment.add(currentPoint);
                                        continue;
                                    }

                                    // 计算当前点与段内最后一个点的距离
                                    TrackPoint lastPoint = currentSegment.get(currentSegment.size() - 1);
                                    double distance = Math.hypot(
                                            currentPoint.getLon() - lastPoint.getLon(),
                                            currentPoint.getLat() - lastPoint.getLat());

                                    if (distance > segmentDistanceThreshold) {
                                        // 距离超过阈值，结束当前段，开始新段
                                        if (currentSegment.size() >= 2) {
                                            segments.add(new ArrayList<>(currentSegment));
                                            log.debug("组{}分段完成: 段内点数={}, 最后距离={}m", k + 1, currentSegment.size(),
                                                    distance);
                                        }
                                        currentSegment.clear();
                                    }

                                    currentSegment.add(currentPoint);
                                }

                                // 处理最后一段
                                if (currentSegment.size() >= 2) {
                                    segments.add(currentSegment);
                                    log.debug("组{}最后一段: 段内点数={}", k + 1, currentSegment.size());
                                }

                                log.debug("组{}轨迹分段完成: 共{}段", k + 1, segments.size());

                                // 对每一段分别构建线缓冲
                                for (int segIndex = 0; segIndex < segments.size(); segIndex++) {
                                    List<TrackPoint> segmentPoints = segments.get(segIndex);

                                    if (segmentPoints.size() < 2) {
                                        log.warn("组{}第{}段点数不足2个，跳过", k + 1, segIndex);
                                        continue;
                                    }

                                    // 构建当前段的线几何
                                    Coordinate[] coords = new Coordinate[segmentPoints.size()];
                                    for (int i = 0; i < segmentPoints.size(); i++) {
                                        TrackPoint pt = segmentPoints.get(i);
                                        coords[i] = new Coordinate(pt.getLon(), pt.getLat());
                                    }

                                    LineString trackLine = config.geometryFactory.createLineString(coords);

                                    // 构建线缓冲（在高斯投影坐标系下计算，单位：米）
                                    Geometry buffer = trackLine.buffer(halfWidth);

                                    log.debug("组{}第{}段线缓冲构建完成: 原始点数={}, 缓冲类型={}, 顶点数={}",
                                            k + 1, segIndex, segmentPoints.size(), buffer.getGeometryType(),
                                            buffer.getNumPoints());

                                    // 处理缓冲结果
                                    if (buffer.isEmpty()) {
                                        log.warn("组{}第{}段线缓冲为空", k + 1, segIndex);
                                    } else if (buffer instanceof Polygon) {
                                        // 单多边形
                                        allPolygons.add(buffer);
                                        log.debug("组{}第{}段添加单多边形", k + 1, segIndex);
                                    } else if (buffer instanceof MultiPolygon) {
                                        // 多多边形
                                        MultiPolygon multiPoly = (MultiPolygon) buffer;
                                        for (int i = 0; i < multiPoly.getNumGeometries(); i++) {
                                            allPolygons.add(multiPoly.getGeometryN(i));
                                            log.debug("组{}第{}段添加子多边形: 索引={}", k + 1, segIndex, i);
                                        }
                                        log.debug("组{}第{}段添加多多边形: 子多边形数量={}", k + 1, segIndex,
                                                multiPoly.getNumGeometries());
                                    } else {
                                        log.warn("组{}第{}段未知的缓冲几何类型: {}", k + 1, segIndex, buffer.getGeometryType());
                                    }
                                }

                                log.debug("组{}所有段处理完成: 总段数={}, 生成子多边形数量={}", k + 1, segments.size(),
                                        allPolygons.size());
                            } catch (Exception e) {
                                log.warn("组{}线缓冲构建失败: {}", k + 1, e.getMessage());
                            }
                        }
                    }
                }
            }
        }

        // 将allPolygons再次进行一次合并
        List<OutlinePart> parts = new ArrayList<>();
        if (CollUtil.isNotEmpty(allPolygons)) {
            try {
                // 使用UnaryUnionOp合并相邻或相交的多边形
                // 先对每个多边形应用缓冲扩展，增加合并的敏感度
                List<Geometry> bufferedPolygons = new ArrayList<>();
                for (Geometry polygon : allPolygons) {
                    // 使用配置中的合并缓冲距离扩展多边形边界
                    Geometry bufferedPolygon = polygon.buffer(halfHalfWidth);
                    bufferedPolygons.add(bufferedPolygon);
                }

                UnaryUnionOp unionOp = new UnaryUnionOp(bufferedPolygons);
                Geometry mergedGeometry = unionOp.union();

                // 合并后再收缩回原始大小（减去缓冲距离）
                if (!mergedGeometry.isEmpty()) {
                    mergedGeometry = mergedGeometry.buffer(-halfHalfWidth);
                }

                // 将合并后的几何转换为几何列表（避免类型强转）
                List<Geometry> mergedGeometries = new ArrayList<>();
                if (mergedGeometry instanceof Polygon || mergedGeometry instanceof MultiPolygon) {
                    if (mergedGeometry instanceof Polygon) {
                        mergedGeometries.add(mergedGeometry);
                    } else {
                        MultiPolygon multiPoly = (MultiPolygon) mergedGeometry;
                        for (int i = 0; i < multiPoly.getNumGeometries(); i++) {
                            mergedGeometries.add(multiPoly.getGeometryN(i));
                        }
                    }
                }

                // 为每个合并后的几何创建OutlinePart
                for (Geometry geom : mergedGeometries) {
                    OutlinePart part = new OutlinePart();
                    part.setTotalWidthM(totalWidthM);
                    part.setOutline(geom);
                    part.setWkt(gaussGeometryToWgs84Wkt(geom, wgs84Points.get(0).getLon()));
                    // 用当前段的所有轨迹点过滤
                    List<TrackPoint> geometryGaussPoints = filterGaussPointsByGeometry(gaussPoints, geom);
                    part.setTrackPoints(
                            gaussPointsToWgs84(geometryGaussPoints, wgs84Points.get(0).getLon()));
                    if (!geometryGaussPoints.isEmpty()) {
                        part.setStartTime(geometryGaussPoints.get(0).getTime());
                        part.setEndTime(
                                geometryGaussPoints.get(geometryGaussPoints.size() - 1).getTime());
                    }
                    part.setMu(calcMuByWgs84Wkt(part.getWkt()));
                    parts.add(part);
                }

                log.info("多边形合并完成: 原始多边形数量={}, 合并后几何数量={}",
                        allPolygons.size(), mergedGeometries.size());

            } catch (Exception e) {
                log.warn("多边形合并失败，使用原始多边形: {}", e.getMessage());
                // 如果合并失败，回退到原始逻辑
                for (Geometry geom : allPolygons) {
                    OutlinePart part = new OutlinePart();
                    part.setTotalWidthM(totalWidthM);
                    part.setOutline(geom);
                    part.setWkt(gaussGeometryToWgs84Wkt(geom, wgs84Points.get(0).getLon()));
                    // 用当前段的所有轨迹点过滤
                    List<TrackPoint> geometryGaussPoints = filterGaussPointsByGeometry(gaussPoints, geom);
                    part.setTrackPoints(
                            gaussPointsToWgs84(geometryGaussPoints, wgs84Points.get(0).getLon()));
                    if (!geometryGaussPoints.isEmpty()) {
                        part.setStartTime(geometryGaussPoints.get(0).getTime());
                        part.setEndTime(
                                geometryGaussPoints.get(geometryGaussPoints.size() - 1).getTime());
                    }
                    part.setMu(calcMuByWgs84Wkt(part.getWkt()));
                    parts.add(part);
                }
            }
        }

        if (CollUtil.isNotEmpty(parts)) {
            double mu = 0.0;
            List<TrackPoint> allTrackPoints = new ArrayList<>();

            // 构建总的MultiPolygon outline
            List<Polygon> polygons = new ArrayList<>();
            for (OutlinePart op : parts) {
                mu += op.getMu();
                allTrackPoints.addAll(op.getTrackPoints());
                polygons.add((Polygon) op.getOutline());
            }

            // 排序所有轨迹点，按时间顺序
            allTrackPoints.sort(Comparator.comparing(TrackPoint::getTime));

            // 创建MultiPolygon
            Geometry finalOutline;
            if (polygons.size() == 1) {
                finalOutline = polygons.get(0);
            } else {
                finalOutline = config.geometryFactory.createMultiPolygon(
                        polygons.toArray(new Polygon[0]));
            }

            result.setParts(parts);
            result.setOutline(finalOutline);
            result.setWkt(gaussGeometryToWgs84Wkt(finalOutline, wgs84Points.get(0).getLon()));
            result.setMu(mu);
            result.setStartTime(allTrackPoints.get(0).getTime());
            result.setEndTime(allTrackPoints.get(allTrackPoints.size() - 1).getTime());
        }

        return result;
    }

    /**
     * WGS84坐标系转换为高斯投影坐标系
     * 
     * @param wgs84Points WGS84坐标系的轨迹点列表
     * @return 高斯投影坐标系的轨迹点列表
     */
    private List<TrackPoint> convertWgs84ToGaussProjection(List<TrackPoint> wgs84Points) {
        if (CollUtil.isEmpty(wgs84Points)) {
            return new ArrayList<>();
        }

        List<TrackPoint> result = new ArrayList<>();

        try {
            // 获取高斯投影坐标系
            CoordinateReferenceSystem gaussCRS = getGaussKrugerCRSByLon(wgs84Points.get(0).getLon());
            CoordinateReferenceSystem wgs84CRS = DefaultGeographicCRS.WGS84;

            // 创建坐标转换
            MathTransform transform = CRS.findMathTransform(wgs84CRS, gaussCRS, true);

            for (TrackPoint point : wgs84Points) {
                try {
                    // 创建WGS84坐标
                    DirectPosition2D wgs84Pos = new DirectPosition2D(wgs84CRS, point.getLon(), point.getLat());

                    // 转换为高斯投影坐标
                    DirectPosition2D gaussPos = new DirectPosition2D();
                    transform.transform(wgs84Pos, gaussPos);

                    // 创建新的轨迹点，使用高斯投影坐标（单位：米）
                    TrackPoint gaussPoint = new TrackPoint(point.getTime(), gaussPos.getX(), gaussPos.getY());

                    result.add(gaussPoint);

                } catch (Exception e) {
                    log.warn("单点坐标转换失败: lat={}, lon={}, error: {}",
                            point.getLat(), point.getLon(), e.getMessage());
                    // 跳过转换失败的点
                }
            }

            log.debug("WGS84转高斯投影完成: 原始点数={}, 转换成功点数={}",
                    wgs84Points.size(), result.size());

        } catch (Exception e) {
            log.warn("坐标转换初始化失败", e.getMessage());
            return new ArrayList<>();
        }

        return result;
    }

    /**
     * 速度过滤：在高斯投影坐标系下计算速度（米单位）
     * 保留相邻点瞬时速度落在开区间 (minSpeedKmh, maxSpeedKmh) 的点
     * 首点始终保留，保持时间升序；边界速度值将被剔除
     *
     * @param gaussPoints 高斯投影坐标系的轨迹点列表（单位：米）
     * @param minSpeedKmh 最小速度阈值（km/h），严格大于该值才保留
     * @param maxSpeedKmh 最大速度阈值（km/h），严格小于该值才保留
     * @return 过滤后的点列表（时间升序）；当输入少于2个点时返回原列表拷贝
     */
    private List<TrackPoint> filterBySpeedRangeMeter(List<TrackPoint> gaussPoints, double minSpeedKmh,
            double maxSpeedKmh) {
        if (gaussPoints == null || gaussPoints.size() <= 1) {
            return gaussPoints == null ? Collections.emptyList() : new ArrayList<>(gaussPoints);
        }
        List<TrackPoint> res = new ArrayList<>(gaussPoints.size());
        res.add(gaussPoints.get(0));
        for (int i = 1; i < gaussPoints.size(); i++) {
            TrackPoint prev = gaussPoints.get(i - 1);
            TrackPoint curr = gaussPoints.get(i);
            double v = speedKmhMeter(prev, curr);
            // 删除小于等于 min 或大于等于 max 的点（删掉 curr）
            if (v > minSpeedKmh && v < maxSpeedKmh) {
                res.add(curr);
            }
        }
        return res;
    }

    /**
     * 在高斯投影坐标系下计算两点间速度（km/h）
     * 使用欧几里得距离（米）与时间差计算瞬时速度
     *
     * @param a 起点（高斯投影坐标，单位：米）
     * @param b 终点（高斯投影坐标，单位：米）
     * @return 速度值（km/h）；无法计算时返回 Double.POSITIVE_INFINITY
     */
    private double speedKmhMeter(TrackPoint a, TrackPoint b) {
        if (a == null || b == null || a.getTime() == null || b.getTime() == null) {
            return Double.POSITIVE_INFINITY;
        }
        // 高斯投影坐标系下使用欧几里得距离（米）
        double dx = b.getLon() - a.getLon(); // X坐标差（东向）
        double dy = b.getLat() - a.getLat(); // Y坐标差（北向）
        double meters = Math.sqrt(dx * dx + dy * dy); // 欧几里得距离（米）

        double seconds = (double) Duration.between(a.getTime(), b.getTime()).toMillis() / 1000.0;
        if (seconds <= 0) {
            return Double.POSITIVE_INFINITY;
        }
        return (meters / seconds) * 3.6; // m/s → km/h
    }

    /**
     * 在高斯投影坐标系下，过滤出在给定几何图形内的轨迹点
     * 渐进式精确过滤：边界框预过滤 → 多线程并行处理 → 100%精确验证
     * 支持多轮廓、复杂几何、确保零误判
     * 
     * @param gaussPoints 高斯投影坐标系的轨迹点列表（单位：米）
     * @param filterGeom  高斯投影坐标系下的过滤几何图形（支持多轮廓）
     * @return 高斯投影坐标系的轨迹点列表（仅包含在filterGeom内的点）
     */
    private List<TrackPoint> filterGaussPointsByGeometry(List<TrackPoint> gaussPoints, Geometry filterGeom) {
        if (CollUtil.isEmpty(gaussPoints) || filterGeom == null || filterGeom.isEmpty()) {
            return new ArrayList<>();
        }

        // ===== 第一阶段：边界框预过滤（内存高效）=====
        Envelope filterEnv = filterGeom.getEnvelopeInternal();
        List<TrackPoint> bboxFiltered = new ArrayList<>();

        for (TrackPoint point : gaussPoints) {
            if (filterEnv.contains(point.getLon(), point.getLat())) {
                bboxFiltered.add(point);
            }
        }

        log.debug("边界框预过滤: 原始点数={}, 边界框内点数={}, 过滤几何类型={}",
                gaussPoints.size(), bboxFiltered.size(), filterGeom.getGeometryType());

        // 如果边界框内无点，直接返回
        if (bboxFiltered.isEmpty()) {
            return bboxFiltered;
        }

        // ===== 第二阶段：多线程精确过滤（确保100%准确）=====
        return preciseParallelFilter(bboxFiltered, filterGeom);
    }

    /**
     * 多线程精确空间过滤（100%准确，零误判）
     * 支持复杂几何、多轮廓、确保每个点都经过精确几何判断
     */
    private List<TrackPoint> preciseParallelFilter(List<TrackPoint> points, Geometry filterGeom) {
        // 小数据集：单线程处理（避免线程开销）
        if (points.size() < 1000) {
            return preciseSingleThreadFilter(points, filterGeom);
        }

        // 大数据集：多线程批处理
        return preciseMultiThreadFilter(points, filterGeom);
    }

    /**
     * 单线程精确过滤（小数据集）
     */
    private List<TrackPoint> preciseSingleThreadFilter(List<TrackPoint> points, Geometry filterGeom) {
        List<TrackPoint> result = new ArrayList<>();

        for (TrackPoint point : points) {
            try {
                Point pointGeom = config.geometryFactory.createPoint(
                        new Coordinate(point.getLon(), point.getLat()));

                // 精确几何判断：支持复杂多轮廓、孔洞等
                if (filterGeom.contains(pointGeom)) {
                    result.add(point);
                }
            } catch (Exception e) {
                log.warn("单线程空间过滤失败: x={}, y={}, error: {}",
                        point.getLon(), point.getLat(), e.getMessage());
            }
        }

        log.trace("单线程精确过滤完成: 输入点数={}, 结果点数={}", points.size(), result.size());
        return result;
    }

    /**
     * 多线程精确过滤（大数据集）
     * 批处理 + 并行流，确保高性能和准确性
     */
    private List<TrackPoint> preciseMultiThreadFilter(List<TrackPoint> points, Geometry filterGeom) {
        final int THREAD_BATCH_SIZE = Math.max(100, points.size() / (Runtime.getRuntime().availableProcessors() * 4));

        List<TrackPoint> result = Collections.synchronizedList(new ArrayList<>());

        // 按批次并行处理
        IntStream.range(0, (points.size() + THREAD_BATCH_SIZE - 1) / THREAD_BATCH_SIZE)
                .parallel()
                .forEach(batchIndex -> {
                    int start = batchIndex * THREAD_BATCH_SIZE;
                    int end = Math.min(start + THREAD_BATCH_SIZE, points.size());
                    List<TrackPoint> batch = points.subList(start, end);

                    // 处理单个批次
                    List<TrackPoint> batchResult = preciseSingleThreadFilter(batch, filterGeom);
                    result.addAll(batchResult);

                    log.trace("批次处理完成: 批次索引={}, 批次大小={}, 结果大小={}",
                            batchIndex, batch.size(), batchResult.size());
                });

        log.debug("多线程精确过滤完成: 输入点数={}, 结果点数={}, 批次大小={}, 处理器数={}",
                points.size(), result.size(), THREAD_BATCH_SIZE, Runtime.getRuntime().availableProcessors());

        return new ArrayList<>(result);
    }

    /**
     * 使用WGS84坐标系的WKT字符串计算面积（亩）
     * 采用与Turf.js相同的球面面积计算算法，确保计算结果一致
     * 
     * @param wkt WGS84坐标系的WKT字符串
     * @return 面积（亩），四舍五入保留4位小数；非法或空输入返回0
     */
    private double calcMuByWgs84Wkt(String wkt) {
        if (StrUtil.isBlank(wkt)) {
            return 0.0;
        }

        try {
            // 解析WKT字符串为几何图形
            WKTReader reader = new WKTReader(config.geometryFactory);
            Geometry geometry = reader.read(wkt);

            if (geometry == null || geometry.isEmpty()) {
                return 0.0;
            }

            // 使用GeoTools的球面面积计算，与Turf.js算法一致
            double areaSqm = calculateSphericalArea(geometry);

            // 转换为亩：1亩 = 2000/3平方米，四舍五入保留4位小数
            return Math.round((areaSqm / (2000.0 / 3.0)) * 10000.0) / 10000.0;

        } catch (Exception e) {
            log.warn("WKT字符串计算面积失败: {}", e.getMessage());
            return 0.0;
        }
    }

    /**
     * 计算球面面积（平方米）
     * 使用GeoTools的球面面积计算算法，与Turf.js的球面面积计算结果一致
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
        } else if (geometry instanceof MultiPolygon) {
            MultiPolygon multiPolygon = (MultiPolygon) geometry;
            for (int i = 0; i < multiPolygon.getNumGeometries(); i++) {
                Polygon polygon = (Polygon) multiPolygon.getGeometryN(i);
                totalArea += calculatePolygonSphericalArea(polygon);
            }
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

        // 减去内环（孔洞）面积
        double holesArea = 0.0;
        for (int i = 0; i < polygon.getNumInteriorRing(); i++) {
            holesArea += calculateRingSphericalArea(polygon.getInteriorRingN(i));
        }

        return exteriorArea - holesArea;
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

        // WGS84椭球长半轴（米），与Turf.js保持一致
        final double EARTH_RADIUS = 6378137.0;

        Coordinate[] coords = ring.getCoordinates();
        if (coords.length < 3) {
            return 0.0; // 需要至少3个点才能形成多边形
        }

        double area = 0.0;

        // 使用球面多边形面积公式
        for (int i = 0; i < coords.length - 1; i++) {
            double lon1 = Math.toRadians(coords[i].x);
            double lat1 = Math.toRadians(coords[i].y);
            double lon2 = Math.toRadians(coords[i + 1].x);
            double lat2 = Math.toRadians(coords[i + 1].y);

            // 球面多边形面积公式
            area += (lon2 - lon1) * (2 + Math.sin(lat1) + Math.sin(lat2));
        }

        area = Math.abs(area) * EARTH_RADIUS * EARTH_RADIUS / 2.0;
        return area;
    }

    /**
     * 计算高斯投影坐标系下几何图形的面积（亩）
     * 直接在高斯投影坐标系（米制）下计算面积，然后转换为亩
     * 1亩 = 666.6667平方米
     * 
     * 算法与 Turf.js 算法返回亩数一致
     * 
     * @param gaussGeometry 高斯投影坐标系下的几何图形（单位：米）
     * @return 面积（亩），四舍五入保留4位小数；非法或空输入返回0
     */
    private double calcMuByGaussGeometry(Geometry gaussGeometry) {
        if (gaussGeometry == null || gaussGeometry.isEmpty()) {
            return 0.0;
        }

        try {
            double areaSqm = 0.0;

            if (gaussGeometry instanceof Polygon) {
                Polygon poly = (Polygon) gaussGeometry;
                // 计算外环面积
                double areaOuter = poly.getArea();
                // 减去内环（孔洞）面积
                double holesArea = 0.0;
                for (int i = 0; i < poly.getNumInteriorRing(); i++) {
                    // 创建内环多边形来计算面积
                    LinearRing hole = poly.getInteriorRingN(i);
                    Polygon holePoly = config.geometryFactory.createPolygon(hole);
                    holesArea += holePoly.getArea();
                }
                areaSqm = areaOuter - holesArea;
            } else if (gaussGeometry instanceof MultiPolygon) {
                MultiPolygon mp = (MultiPolygon) gaussGeometry;
                for (int i = 0; i < mp.getNumGeometries(); i++) {
                    Polygon poly = (Polygon) mp.getGeometryN(i);
                    // 计算每个多边形的外环面积
                    double areaOuter = poly.getArea();
                    // 减去内环（孔洞）面积
                    double holesArea = 0.0;
                    for (int j = 0; j < poly.getNumInteriorRing(); j++) {
                        // 创建内环多边形来计算面积
                        LinearRing hole = poly.getInteriorRingN(j);
                        Polygon holePoly = config.geometryFactory.createPolygon(hole);
                        holesArea += holePoly.getArea();
                    }
                    areaSqm += (areaOuter - holesArea);
                }
            } else {
                return 0.0;
            }

            // 转换为亩：1亩 = 2000/3平方米，四舍五入保留4位小数
            return Math.round((areaSqm / (2000.0 / 3.0)) * 10000.0) / 10000.0;

        } catch (Exception e) {
            log.warn("计算高斯投影面积失败: {}", e.getMessage());
            return 0.0;
        }
    }

    /**
     * 将高斯投影坐标系的轨迹点列表转换为WGS84坐标系
     * 
     * @param gaussPoints 高斯投影坐标系的轨迹点列表（单位：米）
     * @param wgs84Lon    WGS84经度（度），用于确定高斯投影带；通常使用轨迹点的经度
     * @return WGS84坐标系的轨迹点列表；转换失败时返回空列表
     */
    private List<TrackPoint> gaussPointsToWgs84(List<TrackPoint> gaussPoints, double wgs84Lon) {
        if (CollUtil.isEmpty(gaussPoints)) {
            return new ArrayList<>();
        }

        List<TrackPoint> result = new ArrayList<>();

        try {
            // 获取高斯投影坐标系
            CoordinateReferenceSystem gaussCRS = getGaussKrugerCRSByLon(wgs84Lon);
            CoordinateReferenceSystem wgs84CRS = DefaultGeographicCRS.WGS84;

            // 创建坐标转换（高斯投影到WGS84）
            MathTransform transform = CRS.findMathTransform(gaussCRS, wgs84CRS, true);

            for (TrackPoint gaussPoint : gaussPoints) {
                try {
                    // 创建高斯投影坐标
                    DirectPosition2D gaussPos = new DirectPosition2D(gaussCRS, gaussPoint.getLon(),
                            gaussPoint.getLat());

                    // 转换为WGS84坐标
                    DirectPosition2D wgs84Pos = new DirectPosition2D();
                    transform.transform(gaussPos, wgs84Pos);

                    // 创建新的轨迹点，使用WGS84坐标（单位：度）
                    TrackPoint wgs84Point = new TrackPoint(gaussPoint.getTime(), wgs84Pos.getX(), wgs84Pos.getY());

                    result.add(wgs84Point);

                } catch (Exception e) {
                    log.warn("单点坐标转换失败（高斯投影→WGS84）: x={}, y={}, error: {}",
                            gaussPoint.getLon(), gaussPoint.getLat(), e.getMessage());
                    // 跳过转换失败的点
                }
            }

            log.debug("高斯投影转WGS84完成: 原始点数={}, 转换成功点数={}",
                    gaussPoints.size(), result.size());

        } catch (Exception e) {
            log.warn("高斯投影转WGS84坐标转换初始化失败", e);
            return new ArrayList<>();
        }

        return result;
    }

    /**
     * 将高斯投影坐标系的几何图形转换为WGS84坐标系的WKT字符串
     * 
     * @param gaussGeometry 高斯投影坐标系的几何图形
     * @param wgs84Lon      WGS84经度（度），用于确定高斯投影带；通常使用轨迹点的经度
     * @return WGS84坐标系的WKT字符串；转换失败时返回空几何的WKT
     */
    private String gaussGeometryToWgs84Wkt(Geometry gaussGeometry, double wgs84Lon) {
        if (gaussGeometry == null || gaussGeometry.isEmpty()) {
            return config.EMPTY_WKT;
        }

        try {
            // 获取高斯投影到WGS84的坐标转换
            MathTransform txBack = getTxGaussToWgsByLon(wgs84Lon);

            // 将高斯投影几何图形转换到WGS84坐标系
            Geometry wgs84Geometry = JTS.transform(gaussGeometry, txBack);

            // 转换为WKT字符串
            return wgs84Geometry.toText();
        } catch (Exception e) {
            log.warn("高斯投影几何图形转WGS84 WKT失败: {}", e.getMessage());
            return config.EMPTY_WKT;
        }
    }

}