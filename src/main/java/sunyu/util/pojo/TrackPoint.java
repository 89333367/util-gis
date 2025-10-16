package sunyu.util.pojo;

import java.time.LocalDateTime;

/**
 * 轨迹点数据类，继承自CoordinatePoint
 * 包含位置、时间、速度和有效性等信息
 *
 * @author SunYu
 */
public class TrackPoint extends CoordinatePoint {
    /**
     * 时间戳（毫秒）
     */
    private LocalDateTime time;

    /**
     * 速度（km/h）
     */
    private double speed;

    /**
     * 方向角 0~360
     */
    private double direction;


    /**
     * 构造函数，创建一个轨迹点对象
     *
     * @param lon 经度坐标值
     * @param lat 纬度坐标值
     */
    public TrackPoint(double lon, double lat) {
        super(lon, lat);
    }

    public TrackPoint(double lon, double lat, LocalDateTime time, double speed, double direction) {
        super(lon, lat);
        this.time = time;
        this.speed = speed;
        this.direction = direction;
    }

    public LocalDateTime getTime() {
        return time;
    }

    public void setTime(LocalDateTime time) {
        this.time = time;
    }

    public double getSpeed() {
        return speed;
    }

    public void setSpeed(double speed) {
        this.speed = speed;
    }

    public double getDirection() {
        return direction;
    }

    public void setDirection(double direction) {
        this.direction = direction;
    }
}
