package sunyu.util.test;

import sunyu.util.KalmanFilterUtil;
import sunyu.util.pojo.TrackPoint;

import java.time.LocalDateTime;
import java.util.ArrayList;
import java.util.List;

public class KalmanFilterTest {
    public static void main(String[] args) {
        // 1. 创建原始轨迹（仅含经纬度+时间）
        List<TrackPoint> originalTrack = createOriginalTrack();
        System.out.println("==================== 1. 原始轨迹（速度/方向角默认0） ====================");
        printTrack(originalTrack, "原始");

        // 2. 初始化工具类，正确获取预处理后的轨迹（调用工具类暴露的方法）
        KalmanFilterUtil kalmanFilter = new KalmanFilterUtil();
        List<TrackPoint> preprocessedTrack = kalmanFilter.getPreprocessedTrack(originalTrack); // 关键修复：调用public方法
        System.out.println("\n==================== 2. 预处理后轨迹（速度/方向角已计算） ====================");
        printTrack(preprocessedTrack, "预处理后");

        // 3. 执行滤波
        List<TrackPoint> filteredTrack = kalmanFilter.filterTrack(originalTrack);
        System.out.println("\n==================== 3. 滤波后轨迹（异常点已修正） ====================");
        printTrack(filteredTrack, "滤波后");

        // 4. 对比异常点效果
        System.out.println("\n==================== 4. 异常点修正对比 ====================");
        compareAbnormalPoints(preprocessedTrack, filteredTrack);
    }

    /**
     * 创建原始轨迹（含2个异常点）
     */
    private static List<TrackPoint> createOriginalTrack() {
        List<TrackPoint> track = new ArrayList<>();
        LocalDateTime baseTime = LocalDateTime.of(2025, 10, 16, 10, 0, 0);

        for (int i = 0; i < 10; i++) {
            double normalLon = 116.3 + (i * 0.0006); // 每10秒约66米
            double normalLat = 39.9 + (i * 0.0006);

            // 异常点：第4个（i=3）、第8个（i=7）
            if (i == 3) {
                normalLon += 0.002;  // 偏移约222米
                normalLat += 0.0015;
            }
            if (i == 7) {
                normalLon -= 0.0018; // 偏移约199米
                normalLat += 0.0022;
            }

            // 仅传经纬度+时间，速度/方向角默认0.0
            track.add(new TrackPoint(normalLon, normalLat, baseTime.plusSeconds(i * 10)));
        }
        return track;
    }

    /**
     * 打印轨迹详情（正确打印新对象的速度）
     */
    private static void printTrack(List<TrackPoint> track, String type) {
        System.out.printf("%-2s | %-20s | %-12s | %-12s | %-10s | %-10s%n",
                "序号", "时间", "经度(°)", "纬度(°)", "速度(km/h)", "方向角(°)");
        System.out.println("-------------------------------------------------------------------------------------------------");
        for (int i = 0; i < track.size(); i++) {
            TrackPoint p = track.get(i);
            System.out.printf("%-2d | %-20s | %-12.6f | %-12.6f | %-10.1f | %-10.1f%n",
                    i + 1, p.getTime(), p.getLon(), p.getLat(), p.getSpeed(), p.getDirection());
        }
    }

    /**
     * 对比异常点修正效果（用预处理后的轨迹vs滤波后）
     */
    private static void compareAbnormalPoints(List<TrackPoint> preprocessed, List<TrackPoint> filtered) {
        int[] abnormalIndexes = {3, 7}; // 第4个、第8个点（索引3、7）
        for (int idx : abnormalIndexes) {
            TrackPoint pre = preprocessed.get(idx);
            TrackPoint fil = filtered.get(idx);

            double lonOffset = Math.abs(pre.getLon() - fil.getLon());
            double latOffset = Math.abs(pre.getLat() - fil.getLat());

            System.out.printf("第%d个异常点：%n", idx + 1);
            System.out.printf("  预处理后：经纬度(%.6f,%.6f)，速度%.1f km/h%n",
                    pre.getLon(), pre.getLat(), pre.getSpeed());
            System.out.printf("  滤波后：经纬度(%.6f,%.6f)，速度%.1f km/h%n",
                    fil.getLon(), fil.getLat(), fil.getSpeed());
            System.out.printf("  修正偏移：经度%.6f°(%.1f米)，纬度%.6f°(%.1f米)%n%n",
                    lonOffset, lonOffset * 111000, latOffset, latOffset * 111000);
        }
    }
}