package sunyu.util.pojo;

import java.time.LocalDateTime;

/**
 * 轨迹点数据类，继承自CoordinatePoint
 *
 * @author SunYu
 */
public class TrackPoint extends CoordinatePoint {
    /**
     * 定位时间
     */
    private LocalDateTime time;

    public TrackPoint() {
    }

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

    public void setTime(LocalDateTime time) {
        this.time = time;
    }

    @Override
    public double getLat() {
        return super.getLat();
    }

    @Override
    public double getLon() {
        return super.getLon();
    }

    @Override
    public void setLat(double lat) {
        super.setLat(lat);
    }

    @Override
    public void setLon(double lon) {
        super.setLon(lon);
    }

}