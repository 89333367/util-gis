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

    /**
     * GPS速度（单位：m/s）
     */
    private double speed;

    public TrackPoint() {
    }

    /**
     * 完整构造函数：包含所有属性
     * 
     * @param time  时间
     * @param lon   经度
     * @param lat   纬度
     * @param speed GPS速度（单位：m/s）
     */
    public TrackPoint(LocalDateTime time, double lon, double lat, double speed) {
        super(lon, lat);
        this.time = time;
        this.speed = speed;
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

    public double getSpeed() {
        return speed;
    }

    public void setSpeed(double speed) {
        this.speed = speed;
    }

}