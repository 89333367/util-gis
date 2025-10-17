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
     * 时间戳（毫秒）
     */
    private LocalDateTime time;

    /**
     * 速度（km/h）（内部自动计算，用户无需手动设置）
     */
    private double speed;

    /**
     * 方向角 0~360（内部自动计算，用户无需手动设置）
     */
    private double direction;

    // ------------------------------ 新增：用户无需传递speed和direction的构造函数 ------------------------------

    /**
     * 推荐构造函数：用户仅需传递经纬度和时间，速度/方向角由工具类自动计算
     *
     * @param lon  经度
     * @param lat  纬度
     * @param time 时间
     */
    public TrackPoint(double lon, double lat, LocalDateTime time) {
        super(lon, lat);
        this.time = time;
        // 初始值设为0，后续会被工具类覆盖，仅避免null
        this.speed = 0.0;
        this.direction = 0.0;
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

    /**
     * 构造函数：全参数（兼容旧代码，speed/direction会被工具类覆盖）
     *
     * @param lon       经度坐标值
     * @param lat       纬度坐标值
     * @param time      时间
     * @param speed     速度（会被自动计算值覆盖）
     * @param direction 方向角（会被自动计算值覆盖）
     */
    public TrackPoint(double lon, double lat, LocalDateTime time, double speed, double direction) {
        super(lon, lat);
        this.time = time;
        this.speed = speed;
        this.direction = direction;
    }

    // ------------------------------ Getter/Setter（不变） ------------------------------
    public LocalDateTime getTime() {
        return time;
    }

    public void setTime(LocalDateTime time) {
        this.time = time;
    }

    public double getSpeed() {
        return speed;
    }

    // 新增：允许工具类设置计算后的速度
    public void setSpeed(double speed) {
        this.speed = speed;
    }

    public double getDirection() {
        return direction;
    }

    // 新增：允许工具类设置计算后的方向角
    public void setDirection(double direction) {
        this.direction = direction;
    }
}