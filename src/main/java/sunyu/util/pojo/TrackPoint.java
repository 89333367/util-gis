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
     * @param time 时间
     * @param lon  经度
     * @param lat  纬度
     */
    public TrackPoint(LocalDateTime time, double lon, double lat) {
        super(lon, lat);
        this.time = time;
    }

    public LocalDateTime getTime() {
        return time;
    }
    
    public long getTimeMillis() {
        return time.atZone(java.time.ZoneOffset.UTC).toInstant().toEpochMilli();
    }

    public void setTime(LocalDateTime time) {
        this.time = time;
    }

}