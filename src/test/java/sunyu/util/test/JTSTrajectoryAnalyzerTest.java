package sunyu.util.test;


import cn.hutool.log.Log;
import cn.hutool.log.LogFactory;
import org.junit.jupiter.api.Test;
import sunyu.util.JTSTrajectoryAnalyzer;

import java.time.LocalDateTime;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import static org.junit.jupiter.api.Assertions.assertTrue;


/**
 * JTSTrajectoryAnalyzer测试类
 */
public class JTSTrajectoryAnalyzerTest {

    private static final Log logger = LogFactory.get();

    /**
     * 测试梭行法作业轨迹识别
     */
    @Test
    public void testShuttleOperationRecognition() {
        logger.info("开始测试梭行法作业轨迹识别...");

        // 创建分析器
        JTSTrajectoryAnalyzer analyzer = new JTSTrajectoryAnalyzer(10);

        // 生成梭行法作业轨迹
        LocalDateTime baseTime = LocalDateTime.of(2025, 10, 14, 8, 0, 0);
        List<JTSTrajectoryAnalyzer.TrajectoryPoint> shuttleTrajectory = generateShuttleOperation(baseTime, 39.9, 116.3, 100);

        // 分析轨迹
        JTSTrajectoryAnalyzer.AnalysisResult result = analyzer.analyzeTrajectory(shuttleTrajectory);

        // 验证结果
        assertTrue(result.getFieldSegments().size() > 0, "应该识别出田间作业段");
        assertTrue(result.getOverlapDensity() > 1.5, "重叠密度应该大于1.5");
        assertTrue(result.getTotalOperationArea() > 0, "作业面积应该大于0");

        logger.info("梭行法测试通过，识别出{}段田间作业，重叠密度：{}",
                result.getFieldSegments().size(), result.getOverlapDensity());
    }

    /**
     * 测试道路轨迹识别
     */
    @Test
    public void testRoadTrajectoryRecognition() {
        logger.info("开始测试道路轨迹识别...");

        // 创建分析器
        JTSTrajectoryAnalyzer analyzer = new JTSTrajectoryAnalyzer(10);

        // 生成道路轨迹
        LocalDateTime baseTime = LocalDateTime.of(2025, 10, 14, 8, 0, 0);
        List<JTSTrajectoryAnalyzer.TrajectoryPoint> roadTrajectory = generateRoadTrajectory(baseTime, 39.9, 116.3, 50);

        // 分析轨迹
        JTSTrajectoryAnalyzer.AnalysisResult result = analyzer.analyzeTrajectory(roadTrajectory);

        // 验证结果
        assertTrue(result.getRoadSegments().size() > result.getFieldSegments().size(), "应该主要识别为道路行驶");
        assertTrue(result.getOverlapDensity() < 1.2, "重叠密度应该较小");

        logger.info("道路轨迹测试通过，识别出{}段道路行驶", result.getRoadSegments().size());
    }

    /**
     * 测试混合轨迹识别
     */
    @Test
    public void testMixedTrajectoryRecognition() {
        logger.info("开始测试混合轨迹识别...");

        // 创建分析器
        JTSTrajectoryAnalyzer analyzer = new JTSTrajectoryAnalyzer();

        // 生成混合轨迹（道路 + 田间 + 道路）
        LocalDateTime baseTime = LocalDateTime.of(2025, 10, 14, 8, 0, 0);
        List<JTSTrajectoryAnalyzer.TrajectoryPoint> mixedTrajectory = new ArrayList<>();

        // 添加道路段
        mixedTrajectory.addAll(generateRoadTrajectory(baseTime, 39.9, 116.3, 20));
        LocalDateTime nextTime = mixedTrajectory.get(mixedTrajectory.size() - 1).getTime().plusMinutes(2);

        // 添加田间作业段
        mixedTrajectory.addAll(generateShuttleOperation(nextTime,
                mixedTrajectory.get(mixedTrajectory.size() - 1).getLatitude() + 0.01,
                mixedTrajectory.get(mixedTrajectory.size() - 1).getLongitude() + 0.01, 80));
        nextTime = mixedTrajectory.get(mixedTrajectory.size() - 1).getTime().plusMinutes(3);

        // 添加返回道路段
        mixedTrajectory.addAll(generateRoadTrajectory(nextTime,
                mixedTrajectory.get(mixedTrajectory.size() - 1).getLatitude(),
                mixedTrajectory.get(mixedTrajectory.size() - 1).getLongitude(), 20));

        // 分析轨迹
        JTSTrajectoryAnalyzer.AnalysisResult result = analyzer.analyzeTrajectory(mixedTrajectory);

        // 验证结果
        assertTrue(result.getFieldSegments().size() > 0, "应该识别出田间作业段");
        assertTrue(result.getRoadSegments().size() > 0, "应该识别出道路段");

        logger.info("混合轨迹测试通过，识别出{}段田间作业和{}段道路行驶",
                result.getFieldSegments().size(), result.getRoadSegments().size());
    }

    /**
     * 测试性能和内存使用
     */
    @Test
    public void testPerformance() {
        logger.info("开始性能测试...");

        // 创建分析器
        JTSTrajectoryAnalyzer analyzer = new JTSTrajectoryAnalyzer();

        // 生成大规模轨迹数据（10万点）
        LocalDateTime baseTime = LocalDateTime.of(2025, 10, 14, 8, 0, 0);
        List<JTSTrajectoryAnalyzer.TrajectoryPoint> largeTrajectory = generateLargeTrajectory(baseTime, 39.9, 116.3, 100000);

        logger.info("生成了{}个轨迹点", largeTrajectory.size());

        // 测试处理时间
        long startTime = System.currentTimeMillis();
        JTSTrajectoryAnalyzer.AnalysisResult result = analyzer.analyzeTrajectory(largeTrajectory);
        long endTime = System.currentTimeMillis();

        double processingTime = (endTime - startTime) / 1000.0;

        logger.info("处理10万点耗时：{}秒", processingTime);
        assertTrue(processingTime < 60, "处理时间应该小于60秒");

        // 检查内存使用
        Runtime runtime = Runtime.getRuntime();
        long memoryUsed = (runtime.totalMemory() - runtime.freeMemory()) / (1024 * 1024);
        logger.info("内存使用：{}MB", memoryUsed);
        assertTrue(memoryUsed < 200, "内存使用应该小于200MB");
    }

    /**
     * 生成梭行法作业轨迹
     */
    private List<JTSTrajectoryAnalyzer.TrajectoryPoint> generateShuttleOperation(LocalDateTime startTime,
                                                                                 double startLat, double startLon,
                                                                                 int pointCount) {
        List<JTSTrajectoryAnalyzer.TrajectoryPoint> trajectory = new ArrayList<>();
        LocalDateTime currentTime = startTime;
        Random random = new Random(42);

        double currentLat = startLat;
        double currentLon = startLon;
        double mainDirection = 90; // 向东行驶
        boolean isReturning = false;
        double rowSpacing = 0.00009; // 约10米

        for (int i = 0; i < pointCount; i++) {
            double speed = 8 + random.nextDouble() * 4; // 8-12 km/h
            double direction = mainDirection;

            if (isReturning) {
                direction = (direction + 180) % 360;
            }

            // 添加小幅度噪声
            direction += (random.nextDouble() - 0.5) * 2;
            speed += (random.nextDouble() - 0.5) * 0.5;

            // 计算下一个点
            double distance = speed / 3600 * 60; // 1分钟行驶距离
            double latDelta = distance * Math.cos(Math.toRadians(direction)) / 111.32;
            double lonDelta = distance * Math.sin(Math.toRadians(direction)) / 111.32;

            currentLat += latDelta;
            currentLon += lonDelta;

            trajectory.add(new JTSTrajectoryAnalyzer.TrajectoryPoint(
                    currentLat,
                    currentLon,
                    speed,
                    direction,
                    currentTime
            ));

            currentTime = currentTime.plusMinutes(1);

            // 每20个点掉头一次
            if (i > 0 && i % 20 == 0) {
                isReturning = !isReturning;
                // 横向偏移
                double crossDirection = (mainDirection + 90) % 360;
                currentLat += rowSpacing * Math.cos(Math.toRadians(crossDirection));
                currentLon += rowSpacing * Math.sin(Math.toRadians(crossDirection));
            }
        }

        return trajectory;
    }

    /**
     * 生成道路轨迹
     */
    private List<JTSTrajectoryAnalyzer.TrajectoryPoint> generateRoadTrajectory(LocalDateTime startTime,
                                                                               double startLat, double startLon,
                                                                               int pointCount) {
        List<JTSTrajectoryAnalyzer.TrajectoryPoint> trajectory = new ArrayList<>();
        LocalDateTime currentTime = startTime;
        Random random = new Random(43);

        double currentLat = startLat;
        double currentLon = startLon;
        double direction = 45; // 东北方向

        for (int i = 0; i < pointCount; i++) {
            double speed = 40 + random.nextDouble() * 20; // 40-60 km/h

            // 小幅度转向
            if (random.nextDouble() < 0.1) {
                direction += (random.nextDouble() - 0.5) * 5;
                direction = (direction + 360) % 360;
            }

            // 计算下一个点
            double distance = speed / 3600 * 60; // 1分钟行驶距离
            double latDelta = distance * Math.cos(Math.toRadians(direction)) / 111.32;
            double lonDelta = distance * Math.sin(Math.toRadians(direction)) / 111.32;

            currentLat += latDelta;
            currentLon += lonDelta;

            trajectory.add(new JTSTrajectoryAnalyzer.TrajectoryPoint(
                    currentLat,
                    currentLon,
                    speed,
                    direction,
                    currentTime
            ));

            currentTime = currentTime.plusMinutes(1);
        }

        return trajectory;
    }

    /**
     * 生成大规模轨迹数据
     */
    private List<JTSTrajectoryAnalyzer.TrajectoryPoint> generateLargeTrajectory(LocalDateTime startTime,
                                                                                double startLat, double startLon,
                                                                                int pointCount) {
        List<JTSTrajectoryAnalyzer.TrajectoryPoint> trajectory = new ArrayList<>();
        LocalDateTime currentTime = startTime;
        Random random = new Random(44);

        double currentLat = startLat;
        double currentLon = startLon;
        double direction = random.nextDouble() * 360;
        boolean isFieldMode = true;
        int fieldCounter = 0;

        for (int i = 0; i < pointCount; i++) {
            double speed;
            double directionChange;

            if (isFieldMode) {
                // 田间作业模式
                speed = 6 + random.nextDouble() * 6; // 6-12 km/h
                directionChange = (random.nextDouble() - 0.5) * 3;

                // 定期掉头
                if (i % 500 == 0) {
                    direction = (direction + 180) % 360;
                }

                fieldCounter++;
                if (fieldCounter > 10000) {
                    isFieldMode = false;
                    fieldCounter = 0;
                }
            } else {
                // 道路模式
                speed = 30 + random.nextDouble() * 30; // 30-60 km/h
                directionChange = (random.nextDouble() - 0.5) * 2;

                fieldCounter++;
                if (fieldCounter > 2000) {
                    isFieldMode = true;
                    fieldCounter = 0;
                }
            }

            direction += directionChange;
            direction = (direction + 360) % 360;

            // 计算下一个点
            double distance = speed / 3600 * 30; // 30秒行驶距离
            double latDelta = distance * Math.cos(Math.toRadians(direction)) / 111.32;
            double lonDelta = distance * Math.sin(Math.toRadians(direction)) / 111.32;

            currentLat += latDelta;
            currentLon += lonDelta;

            trajectory.add(new JTSTrajectoryAnalyzer.TrajectoryPoint(
                    currentLat,
                    currentLon,
                    speed,
                    direction,
                    currentTime
            ));

            currentTime = currentTime.plusSeconds(30);
        }

        return trajectory;
    }
}
