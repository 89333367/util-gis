package sunyu.util;

import cn.hutool.log.Log;
import cn.hutool.log.LogFactory;
import org.locationtech.jts.geom.*;
import org.locationtech.jts.operation.buffer.BufferOp;
import org.locationtech.jts.operation.buffer.BufferParameters;
import org.locationtech.jts.operation.union.UnaryUnionOp;

import java.time.LocalDateTime;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.Objects;
import java.util.stream.Collectors;

/**
 * 基于JTS的田间作业轨迹分析器
 * 核心算法：通过缓冲区重叠分析识别田间作业轨迹
 */
public class JTSTrajectoryAnalyzer {

    private static final Log logger = LogFactory.get();

    // 默认配置参数
    private static final double DEFAULT_WORKING_WIDTH = 2.5;        // 默认作业幅宽 米
    private static final double DEFAULT_OVERLAP_THRESHOLD = 0.3;     // 重叠阈值 30%
    private static final double DEFAULT_MIN_SEGMENT_LENGTH = 5.0;    // 最小轨迹段长度 米
    private static final double DEFAULT_SPEED_THRESHOLD = 15.0;      // 速度阈值 km/h
    private static final int DEFAULT_MIN_POINTS_PER_SEGMENT = 5;     // 每段最小点数
    private static final double DEFAULT_SPATIAL_TOLERANCE = 0.00001; // 空间容差 度

    // JTS几何工厂
    private final GeometryFactory geometryFactory;

    // 配置参数
    private double workingWidth;
    private double overlapThreshold;
    private double minSegmentLength;
    private double speedThreshold;
    private int minPointsPerSegment;
    private double spatialTolerance;

    /**
     * 轨迹点数据结构
     */
    public static class TrajectoryPoint {
        private double latitude;      // 纬度 (度)
        private double longitude;     // 经度 (度)
        private double speed;         // 速度 (km/h)
        private double direction;     // 方向角 (度，0-360)
        private LocalDateTime time;   // 时间戳

        private boolean isFieldPoint; // 是否为田间作业点
        private double distanceToNext; // 到下一个点的距离 (米)

        public TrajectoryPoint(double latitude, double longitude, double speed,
                               double direction, LocalDateTime time) {
            this.latitude = latitude;
            this.longitude = longitude;
            this.speed = speed;
            this.direction = direction;
            this.time = time;
        }

        // getters and setters
        public double getLatitude() {
            return latitude;
        }

        public double getLongitude() {
            return longitude;
        }

        public double getSpeed() {
            return speed;
        }

        public double getDirection() {
            return direction;
        }

        public LocalDateTime getTime() {
            return time;
        }

        public boolean isFieldPoint() {
            return isFieldPoint;
        }

        public void setFieldPoint(boolean fieldPoint) {
            isFieldPoint = fieldPoint;
        }

        public double getDistanceToNext() {
            return distanceToNext;
        }

        public void setDistanceToNext(double distanceToNext) {
            this.distanceToNext = distanceToNext;
        }

        @Override
        public String toString() {
            return String.format("Point{lat=%.6f, lon=%.6f, speed=%.1f, time=%s, isField=%b}",
                    latitude, longitude, speed, time, isFieldPoint);
        }
    }

    /**
     * 轨迹段数据结构
     */
    public static class TrajectorySegment {
        private List<TrajectoryPoint> points;
        private SegmentType type;
        private LineString lineString;
        private Polygon bufferPolygon;
        private double segmentLength;
        private double overlapRatio;
        private double averageSpeed;

        public enum SegmentType {
            FIELD_OPERATION,  // 田间作业
            ROAD_TRAVEL,      // 道路行驶
            UNKNOWN           // 未知
        }

        public TrajectorySegment(List<TrajectoryPoint> points) {
            this.points = new ArrayList<>(points);
            this.type = SegmentType.UNKNOWN;
            calculateSegmentStats();
        }

        private void calculateSegmentStats() {
            if (points.isEmpty()) return;

            // 计算平均速度
            averageSpeed = points.stream()
                    .mapToDouble(TrajectoryPoint::getSpeed)
                    .average()
                    .orElse(0);

            // 计算总长度
            segmentLength = points.stream()
                    .mapToDouble(TrajectoryPoint::getDistanceToNext)
                    .sum();
        }

        // getters and setters
        public List<TrajectoryPoint> getPoints() {
            return new ArrayList<>(points);
        }

        public SegmentType getType() {
            return type;
        }

        public void setType(SegmentType type) {
            this.type = type;
        }

        public LineString getLineString() {
            return lineString;
        }

        public void setLineString(LineString lineString) {
            this.lineString = lineString;
        }

        public Polygon getBufferPolygon() {
            return bufferPolygon;
        }

        public void setBufferPolygon(Polygon bufferPolygon) {
            this.bufferPolygon = bufferPolygon;
        }

        public double getSegmentLength() {
            return segmentLength;
        }

        public double getOverlapRatio() {
            return overlapRatio;
        }

        public void setOverlapRatio(double overlapRatio) {
            this.overlapRatio = overlapRatio;
        }

        public double getAverageSpeed() {
            return averageSpeed;
        }

        @Override
        public String toString() {
            return String.format("Segment{type=%s, points=%d, length=%.1fm, overlap=%.1f%%, avgSpeed=%.1fkm/h}",
                    type, points.size(), segmentLength, overlapRatio * 100, averageSpeed);
        }
    }

    /**
     * 分析结果数据结构
     */
    public static class AnalysisResult {
        private List<TrajectorySegment> allSegments;
        private List<TrajectorySegment> fieldSegments;
        private List<TrajectorySegment> roadSegments;
        private Polygon totalFieldBuffer;
        private double totalFieldArea;
        private double totalOperationArea;
        private double overlapDensity;
        private double workingWidth;

        public AnalysisResult(List<TrajectorySegment> allSegments, double workingWidth) {
            this.allSegments = allSegments;
            this.workingWidth = workingWidth;
            this.fieldSegments = allSegments.stream()
                    .filter(seg -> seg.getType() == TrajectorySegment.SegmentType.FIELD_OPERATION)
                    .collect(Collectors.toList());
            this.roadSegments = allSegments.stream()
                    .filter(seg -> seg.getType() == TrajectorySegment.SegmentType.ROAD_TRAVEL)
                    .collect(Collectors.toList());

            calculateStats();
        }

        private void calculateStats() {
            // 计算总田间缓冲区
            List<Geometry> fieldBuffers = fieldSegments.stream()
                    .map(TrajectorySegment::getBufferPolygon)
                    .filter(Objects::nonNull)
                    .collect(Collectors.toList());

            if (!fieldBuffers.isEmpty()) {
                totalFieldBuffer = (Polygon) UnaryUnionOp.union(fieldBuffers);
                totalFieldArea = totalFieldBuffer.getArea() * 100000000; // 转换为平方米 (1度约=111319.9米)
            }

            // 计算总作业面积（基于轨迹长度和幅宽）
            totalOperationArea = fieldSegments.stream()
                    .mapToDouble(seg -> seg.getSegmentLength() * workingWidth / 10000) // 转换为公顷
                    .sum();

            // 计算重叠密度
            if (totalFieldArea > 0) {
                double totalBufferArea = fieldSegments.stream()
                        .mapToDouble(seg -> seg.getBufferPolygon().getArea() * 100000000)
                        .sum();
                overlapDensity = totalBufferArea / totalFieldArea;
            }
        }

        // getters
        public List<TrajectorySegment> getAllSegments() {
            return allSegments;
        }

        public List<TrajectorySegment> getFieldSegments() {
            return fieldSegments;
        }

        public List<TrajectorySegment> getRoadSegments() {
            return roadSegments;
        }

        public Polygon getTotalFieldBuffer() {
            return totalFieldBuffer;
        }

        public double getTotalFieldArea() {
            return totalFieldArea;
        }

        public double getTotalOperationArea() {
            return totalOperationArea;
        }

        public double getOverlapDensity() {
            return overlapDensity;
        }

        @Override
        public String toString() {
            return String.format("AnalysisResult{fieldSegments=%d, roadSegments=%d, totalArea=%.2f公顷, overlapDensity=%.2f}",
                    fieldSegments.size(), roadSegments.size(), totalOperationArea, overlapDensity);
        }
    }

    /**
     * 构造函数
     */
    public JTSTrajectoryAnalyzer() {
        this.geometryFactory = new GeometryFactory(new PrecisionModel(PrecisionModel.FLOATING), 4326); // WGS84坐标系
        this.workingWidth = DEFAULT_WORKING_WIDTH;
        this.overlapThreshold = DEFAULT_OVERLAP_THRESHOLD;
        this.minSegmentLength = DEFAULT_MIN_SEGMENT_LENGTH;
        this.speedThreshold = DEFAULT_SPEED_THRESHOLD;
        this.minPointsPerSegment = DEFAULT_MIN_POINTS_PER_SEGMENT;
        this.spatialTolerance = DEFAULT_SPATIAL_TOLERANCE;
    }

    public JTSTrajectoryAnalyzer(double workingWidth) {
        this.geometryFactory = new GeometryFactory(new PrecisionModel(PrecisionModel.FLOATING), 4326); // WGS84坐标系
        this.workingWidth = workingWidth;
        this.overlapThreshold = DEFAULT_OVERLAP_THRESHOLD;
        this.minSegmentLength = DEFAULT_MIN_SEGMENT_LENGTH;
        this.speedThreshold = DEFAULT_SPEED_THRESHOLD;
        this.minPointsPerSegment = DEFAULT_MIN_POINTS_PER_SEGMENT;
        this.spatialTolerance = DEFAULT_SPATIAL_TOLERANCE;
    }

    /**
     * 构造函数 - 自定义参数
     */
    public JTSTrajectoryAnalyzer(double workingWidth, double overlapThreshold,
                                 double minSegmentLength, double speedThreshold,
                                 int minPointsPerSegment, double spatialTolerance) {
        this.geometryFactory = new GeometryFactory(new PrecisionModel(PrecisionModel.FLOATING), 4326);
        this.workingWidth = workingWidth;
        this.overlapThreshold = overlapThreshold;
        this.minSegmentLength = minSegmentLength;
        this.speedThreshold = speedThreshold;
        this.minPointsPerSegment = minPointsPerSegment;
        this.spatialTolerance = spatialTolerance;
    }

    /**
     * 分析轨迹数据
     *
     * @param trajectoryPoints 原始轨迹点列表
     *
     * @return 分析结果
     */
    public AnalysisResult analyzeTrajectory(List<TrajectoryPoint> trajectoryPoints) {
        if (trajectoryPoints == null || trajectoryPoints.size() < minPointsPerSegment) {
            throw new IllegalArgumentException("轨迹点数量不足，至少需要" + minPointsPerSegment + "个点");
        }

        long startTime = System.currentTimeMillis();
        logger.info("开始轨迹分析，总点数：{}", trajectoryPoints.size());

        try {
            // 1. 数据预处理
            List<TrajectoryPoint> processedPoints = preprocessData(trajectoryPoints);
            logger.info("数据预处理完成，剩余点数：{}", processedPoints.size());

            // 2. 轨迹分段
            List<TrajectorySegment> segments = segmentTrajectory(processedPoints);
            logger.info("轨迹分段完成，段数：{}", segments.size());

            // 3. 生成几何对象
            createGeometricObjects(segments);
            logger.info("几何对象生成完成");

            // 4. 缓冲区分析
            generateBuffers(segments);
            logger.info("缓冲区生成完成");

            // 5. 重叠分析
            analyzeOverlaps(segments);
            logger.info("重叠分析完成");

            // 6. 轨迹分类
            classifySegments(segments);
            logger.info("轨迹分类完成");

            // 7. 生成分析结果
            AnalysisResult result = new AnalysisResult(segments, this.workingWidth);
            logger.info("分析完成，耗时：{}ms", System.currentTimeMillis() - startTime);

            return result;

        } catch (Exception e) {
            logger.error("轨迹分析失败", e);
            throw new RuntimeException("轨迹分析失败", e);
        }
    }

    /**
     * 数据预处理
     */
    private List<TrajectoryPoint> preprocessData(List<TrajectoryPoint> rawPoints) {
        // 1. 按时间排序
        List<TrajectoryPoint> sortedPoints = new ArrayList<>(rawPoints);
        sortedPoints.sort(Comparator.comparing(TrajectoryPoint::getTime));

        // 2. 去除重复点和异常点
        List<TrajectoryPoint> cleanedPoints = new ArrayList<>();
        TrajectoryPoint previous = null;

        for (TrajectoryPoint current : sortedPoints) {
            if (previous == null) {
                cleanedPoints.add(current);
                previous = current;
                continue;
            }

            // 检查时间顺序
            if (current.getTime().isBefore(previous.getTime())) {
                logger.warn("时间顺序异常，跳过点：{}", current);
                continue;
            }

            // 检查距离是否过近（重复点）
            double distance = calculateDistance(previous, current);
            if (distance < spatialTolerance * 111319.9) { // 转换为米
                continue;
            }

            // 检查速度是否异常
            if (current.getSpeed() < 0 || current.getSpeed() > 120) {
                logger.warn("速度异常，跳过点：{}", current);
                continue;
            }

            // 计算到下一个点的距离
            previous.setDistanceToNext(distance);

            cleanedPoints.add(current);
            previous = current;
        }

        return cleanedPoints;
    }

    /**
     * 轨迹分段
     */
    private List<TrajectorySegment> segmentTrajectory(List<TrajectoryPoint> points) {
        List<TrajectorySegment> segments = new ArrayList<>();
        if (points.size() < minPointsPerSegment) return segments;

        List<TrajectoryPoint> currentSegmentPoints = new ArrayList<>();
        currentSegmentPoints.add(points.get(0));

        double cumulativeLength = 0;
        TrajectoryPoint previous = points.get(0);

        for (int i = 1; i < points.size(); i++) {
            TrajectoryPoint current = points.get(i);
            currentSegmentPoints.add(current);

            // 累加长度
            cumulativeLength += current.getDistanceToNext();

            // 检查是否需要分段
            boolean needSegmentation = false;

            // 长度达到阈值
            if (cumulativeLength >= minSegmentLength && currentSegmentPoints.size() >= minPointsPerSegment) {
                needSegmentation = true;
            }

            // 速度变化较大
            double speedChange = Math.abs(current.getSpeed() - previous.getSpeed());
            if (speedChange > Math.max(5, previous.getSpeed() * 0.5)) {
                needSegmentation = true;
            }

            // 方向变化较大
            double directionChange = Math.abs(current.getDirection() - previous.getDirection());
            directionChange = Math.min(directionChange, 360 - directionChange);
            if (directionChange > 45) {
                needSegmentation = true;
            }

            if (needSegmentation) {
                segments.add(new TrajectorySegment(currentSegmentPoints));
                currentSegmentPoints.clear();
                currentSegmentPoints.add(current);
                cumulativeLength = 0;
            }

            previous = current;
        }

        // 添加最后一段
        if (currentSegmentPoints.size() >= minPointsPerSegment && cumulativeLength >= minSegmentLength * 0.5) {
            segments.add(new TrajectorySegment(currentSegmentPoints));
        }

        return segments;
    }

    /**
     * 创建几何对象
     */
    private void createGeometricObjects(List<TrajectorySegment> segments) {
        for (TrajectorySegment segment : segments) {
            List<TrajectoryPoint> points = segment.getPoints();
            if (points.size() < 2) continue;

            // 创建坐标序列
            Coordinate[] coordinates = new Coordinate[points.size()];
            for (int i = 0; i < points.size(); i++) {
                TrajectoryPoint point = points.get(i);
                coordinates[i] = new Coordinate(point.getLongitude(), point.getLatitude());
            }

            // 创建线串
            LineString lineString = geometryFactory.createLineString(coordinates);
            segment.setLineString(lineString);
        }
    }

    /**
     * 生成缓冲区
     */
    private void generateBuffers(List<TrajectorySegment> segments) {
        // 将作业幅宽转换为度（假设在赤道附近，1度约=111319.9米）
        double bufferDistanceInDegrees = workingWidth / 2 / 111319.9;

        // 配置缓冲区参数
        BufferParameters bufferParams = new BufferParameters();
        bufferParams.setEndCapStyle(BufferParameters.CAP_ROUND);
        bufferParams.setJoinStyle(BufferParameters.JOIN_ROUND);
        bufferParams.setQuadrantSegments(8); // 圆弧精度

        for (TrajectorySegment segment : segments) {
            LineString lineString = segment.getLineString();
            if (lineString == null) continue;

            try {
                // 生成缓冲区
                Geometry buffer = BufferOp.bufferOp(lineString, bufferDistanceInDegrees, bufferParams);

                if (buffer instanceof Polygon) {
                    segment.setBufferPolygon((Polygon) buffer);
                } else if (buffer instanceof MultiPolygon) {
                    // 如果是多个多边形，取面积最大的一个
                    MultiPolygon multiPolygon = (MultiPolygon) buffer;
                    Polygon largestPolygon = null;
                    double maxArea = 0;

                    for (int i = 0; i < multiPolygon.getNumGeometries(); i++) {
                        Polygon polygon = (Polygon) multiPolygon.getGeometryN(i);
                        if (polygon.getArea() > maxArea) {
                            maxArea = polygon.getArea();
                            largestPolygon = polygon;
                        }
                    }

                    segment.setBufferPolygon(largestPolygon);
                }

            } catch (Exception e) {
                logger.error("生成缓冲区失败", e);
            }
        }
    }

    /**
     * 重叠分析
     */
    private void analyzeOverlaps(List<TrajectorySegment> segments) {
        // 首先合并所有缓冲区用于重叠分析
        List<Geometry> allBuffers = segments.stream()
                .map(TrajectorySegment::getBufferPolygon)
                .filter(Objects::nonNull)
                .collect(Collectors.toList());

        if (allBuffers.isEmpty()) return;

        // 计算总合并缓冲区
        Geometry totalUnion = UnaryUnionOp.union(allBuffers);
        double totalUnionArea = totalUnion.getArea();

        // 计算每个段的重叠比例
        for (TrajectorySegment segment : segments) {
            Polygon buffer = segment.getBufferPolygon();
            if (buffer == null) continue;

            try {
                // 计算与总合并缓冲区的重叠
                Geometry intersection = buffer.intersection(totalUnion);
                double overlapArea = intersection.getArea();
                double bufferArea = buffer.getArea();

                // 重叠比例 = 重叠面积 / 缓冲区面积
                double overlapRatio = bufferArea > 0 ? overlapArea / bufferArea : 0;
                segment.setOverlapRatio(overlapRatio);

            } catch (Exception e) {
                logger.error("计算重叠失败", e);
                segment.setOverlapRatio(0);
            }
        }
    }

    /**
     * 轨迹分类
     */
    private void classifySegments(List<TrajectorySegment> segments) {
        for (TrajectorySegment segment : segments) {
            // 基于多个特征进行分类
            double fieldScore = calculateFieldScore(segment);

            if (fieldScore >= 0.6) {
                segment.setType(TrajectorySegment.SegmentType.FIELD_OPERATION);
                // 标记段内所有点为田间作业点
                segment.getPoints().forEach(point -> point.setFieldPoint(true));
            } else {
                segment.setType(TrajectorySegment.SegmentType.ROAD_TRAVEL);
                // 标记段内所有点为道路点
                segment.getPoints().forEach(point -> point.setFieldPoint(false));
            }
        }
    }

    /**
     * 计算田间作业得分
     */
    private double calculateFieldScore(TrajectorySegment segment) {
        double score = 0;

        // 1. 重叠比例特征（最重要）
        if (segment.getOverlapRatio() >= overlapThreshold) {
            score += 0.4; // 重叠比例高，权重最大
        } else if (segment.getOverlapRatio() > 0) {
            score += segment.getOverlapRatio() * 0.4 / overlapThreshold;
        }

        // 2. 速度特征
        if (segment.getAverageSpeed() >= 1 && segment.getAverageSpeed() <= speedThreshold) {
            score += 0.2;
        }

        // 3. 长度特征
        if (segment.getSegmentLength() > minSegmentLength * 2) {
            score += 0.1;
        }

        // 4. 缓冲区质量特征
        if (segment.getBufferPolygon() != null) {
            score += 0.1;
        }

        // 5. 点数特征
        if (segment.getPoints().size() > minPointsPerSegment * 2) {
            score += 0.1;
        }

        // 6. 速度稳定性特征
        double speedStd = calculateSpeedStandardDeviation(segment.getPoints());
        if (speedStd / segment.getAverageSpeed() < 0.3 && segment.getAverageSpeed() > 0) {
            score += 0.1;
        }

        return score;
    }

    /**
     * 计算速度标准差
     */
    private double calculateSpeedStandardDeviation(List<TrajectoryPoint> points) {
        List<Double> speeds = points.stream()
                .mapToDouble(TrajectoryPoint::getSpeed)
                .boxed()
                .collect(Collectors.toList());

        double mean = speeds.stream().mapToDouble(Double::doubleValue).average().orElse(0);
        double variance = speeds.stream()
                .mapToDouble(s -> Math.pow(s - mean, 2))
                .average()
                .orElse(0);

        return Math.sqrt(variance);
    }

    /**
     * 计算两点之间的距离（米）
     */
    private double calculateDistance(TrajectoryPoint p1, TrajectoryPoint p2) {
        double lat1 = Math.toRadians(p1.getLatitude());
        double lon1 = Math.toRadians(p1.getLongitude());
        double lat2 = Math.toRadians(p2.getLatitude());
        double lon2 = Math.toRadians(p2.getLongitude());

        double dLat = lat2 - lat1;
        double dLon = lon2 - lon1;

        double a = Math.sin(dLat / 2) * Math.sin(dLat / 2) +
                Math.cos(lat1) * Math.cos(lat2) *
                        Math.sin(dLon / 2) * Math.sin(dLon / 2);
        double c = 2 * Math.atan2(Math.sqrt(a), Math.sqrt(1 - a));

        return 6371000 * c; // 地球半径6371000米
    }

    // getters and setters
    public double getWorkingWidth() {
        return workingWidth;
    }

    public void setWorkingWidth(double workingWidth) {
        this.workingWidth = workingWidth;
    }

    public double getOverlapThreshold() {
        return overlapThreshold;
    }

    public void setOverlapThreshold(double overlapThreshold) {
        this.overlapThreshold = overlapThreshold;
    }

    public double getMinSegmentLength() {
        return minSegmentLength;
    }

    public void setMinSegmentLength(double minSegmentLength) {
        this.minSegmentLength = minSegmentLength;
    }

    public double getSpeedThreshold() {
        return speedThreshold;
    }

    public void setSpeedThreshold(double speedThreshold) {
        this.speedThreshold = speedThreshold;
    }

    /**
     * 工具方法：将轨迹转换为WKT格式
     */
    public String convertToWKT(List<TrajectoryPoint> points) {
        if (points.size() < 2) return "";

        StringBuilder wkt = new StringBuilder("LINESTRING (");
        for (int i = 0; i < points.size(); i++) {
            TrajectoryPoint point = points.get(i);
            if (i > 0) wkt.append(", ");
            wkt.append(point.getLongitude()).append(" ").append(point.getLatitude());
        }
        wkt.append(")");

        return wkt.toString();
    }

    /**
     * 工具方法：计算作业面积
     */
    public double calculateOperationArea(List<TrajectorySegment> fieldSegments) {
        return fieldSegments.stream()
                .mapToDouble(seg -> seg.getSegmentLength() * workingWidth / 10000) // 转换为公顷
                .sum();
    }
}
