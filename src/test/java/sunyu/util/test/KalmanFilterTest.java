package sunyu.util.test;

import sunyu.util.KalmanFilterUtil;
import sunyu.util.pojo.TrackPoint;

import java.time.LocalDateTime;
import java.util.ArrayList;
import java.util.List;

/**
 * 卡尔曼滤波工具使用测试类
 */
public class KalmanFilterTest {
    public static void main(String[] args) {
        // 1. 第一步：创建模拟轨迹数据（包含正常点+异常点，模拟真实定位场景）
        List<TrackPoint> originalTrack = createSimulateTrackData();
        System.out.println("=== 原始轨迹点数量：" + originalTrack.size() + " ===");
        printTrackPoints("原始轨迹", originalTrack);

        // 2. 第二步：初始化卡尔曼滤波工具
        KalmanFilterUtil kalmanFilter = new KalmanFilterUtil();

        // 3. 第三步：执行轨迹滤波（核心调用）
        List<TrackPoint> filteredTrack = kalmanFilter.filterTrack(originalTrack);
        System.out.println("\n=== 滤波后轨迹点数量：" + filteredTrack.size() + " ===");
        printTrackPoints("滤波后轨迹", filteredTrack);

        // 4. 第四步：对比关键异常点的滤波效果
        compareAbnormalPoints(originalTrack, filteredTrack);
    }

    /**
     * 模拟创建轨迹数据：沿「东经116.3°→116.4°，北纬39.9°→40.0°」移动，中间插入2个异常点
     *
     * @return 包含正常点和异常点的轨迹列表
     */
    private static List<TrackPoint> createSimulateTrackData() {
        List<TrackPoint> trackPoints = new ArrayList<>();
        LocalDateTime baseTime = LocalDateTime.of(2025, 10, 16, 10, 0, 0); // 基础时间

        // 生成10个连续点，其中第4个、第8个是异常点（故意偏移经纬度）
        for (int i = 0; i < 10; i++) {
            // 正常点规律：每10秒移动0.01°经度、0.01°纬度
            double normalLon = 116.3 + (i * 0.01); // 经度：116.3 → 116.4
            double normalLat = 39.9 + (i * 0.01); // 纬度：39.9 → 40.0

            // 速度：模拟车辆行驶，30~40 km/h
            double speed = 30 + (Math.random() * 10);
            // 方向角：大致东北方向（45°左右，微小波动）
            double direction = 45 + (Math.random() * 5) - 2.5;

            // 插入异常点（第4个点：i=3；第8个点：i=7）
            double actualLon = normalLon;
            double actualLat = normalLat;
            if (i == 3) {
                actualLon += 0.05; // 经度异常偏移+0.05°（约5公里）
                actualLat += 0.04; // 纬度异常偏移+0.04°
            }
            if (i == 7) {
                actualLon -= 0.03; // 经度异常偏移-0.03°
                actualLat += 0.06; // 纬度异常偏移+0.06°
            }

            // 添加到轨迹列表（时间每步加10秒）
            trackPoints.add(new TrackPoint(
                    actualLon,
                    actualLat,
                    baseTime.plusSeconds(i * 10), // 时间递增
                    speed,
                    direction
            ));
        }
        return trackPoints;
    }

    /**
     * 打印轨迹点详情（经纬度、时间、速度）
     *
     * @param trackName   轨迹名称（用于区分原始/滤波后）
     * @param trackPoints 轨迹点列表
     */
    private static void printTrackPoints(String trackName, List<TrackPoint> trackPoints) {
        for (int i = 0; i < trackPoints.size(); i++) {
            TrackPoint point = trackPoints.get(i);
            System.out.printf(
                    "%s 第%d点：经度=%.6f，纬度=%.6f，时间=%s，速度=%.1f km/h\n",
                    trackName,
                    (i + 1),
                    point.getLon(),
                    point.getLat(),
                    point.getTime(),
                    point.getSpeed()
            );
        }
    }

    /**
     * 对比原始轨迹与滤波轨迹的异常点修正效果
     *
     * @param original 原始轨迹
     * @param filtered 滤波后轨迹
     */
    private static void compareAbnormalPoints(List<TrackPoint> original, List<TrackPoint> filtered) {
        System.out.println("\n=== 异常点滤波效果对比 ===");
        // 只对比之前插入的2个异常点（索引3和7，对应第4、第8个点）
        int[] abnormalIndexes = {3, 7};
        for (int index : abnormalIndexes) {
            TrackPoint originalPoint = original.get(index);
            TrackPoint filteredPoint = filtered.get(index);

            // 计算经纬度偏移修正量（单位：度，1度≈111公里）
            double lonCorrection = Math.abs(originalPoint.getLon() - filteredPoint.getLon());
            double latCorrection = Math.abs(originalPoint.getLat() - filteredPoint.getLat());

            System.out.printf(
                    "第%d个异常点：\n" +
                            "  原始经纬度：(%.6f, %.6f)\n" +
                            "  滤波后经纬度：(%.6f, %.6f)\n" +
                            "  修正偏移：经度%.6f°，纬度%.6f°（约合经度%.1f米，纬度%.1f米）\n",
                    (index + 1),
                    originalPoint.getLon(), originalPoint.getLat(),
                    filteredPoint.getLon(), filteredPoint.getLat(),
                    lonCorrection, latCorrection,
                    lonCorrection * 111000, // 1度≈111000米
                    latCorrection * 111000
            );
        }
    }
}