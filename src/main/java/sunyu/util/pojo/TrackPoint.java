package sunyu.util.pojo;

import java.time.LocalDateTime;

/**
 * 轨迹点数据类，继承自CoordinatePoint
 * 包含位置、时间、速度和有效性等信息（速度/方向角可自动计算，无需用户传递）
 *
 * @author SunYu
 */
public class TrackPoint extends CoordinatePoint {
    /**
     * 定位时间
     */
    private LocalDateTime time;

    /**
     * 推荐构造函数：用户仅需传递经纬度和时间
     *
     * @param lon  经度
     * @param lat  纬度
     * @param time 时间
     */
    public TrackPoint(double lon, double lat, LocalDateTime time) {
        super(lon, lat);
        this.time = time;
    }

    // ------------------------------ 保留原有构造函数（保证兼容性） ------------------------------

    /**
     * 构造函数：仅传经纬度（不推荐，时间为null会被工具类过滤）
     *
     * @param lon 经度坐标值
     * @param lat 纬度坐标值
     */
    public TrackPoint(double lon, double lat) {
        super(lon, lat);
    }

    // ------------------------------ Getter/Setter（不变） ------------------------------
    public LocalDateTime getTime() {
        return time;
    }

    public void setTime(LocalDateTime time) {
        this.time = time;
    }

}