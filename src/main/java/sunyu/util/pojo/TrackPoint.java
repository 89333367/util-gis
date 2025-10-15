package sunyu.util.pojo;

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
    private long ts;

    /**
     * 速度（km/h）
     */
    private double speed;

    /**
     * 方向角 0~360
     */
    private int heading;

    /**
     * 有效性标记
     */
    private boolean valid;

    /**
     * 构造函数
     *
     * @param lon   经度
     * @param lat   纬度
     * @param ts    时间戳
     * @param speed 速度
     * @param valid 有效性标记
     */
    public TrackPoint(double lon, double lat, long ts, double speed, boolean valid) {
        super(lon, lat);
        this.ts = ts;
        this.speed = speed;
        this.valid = valid;
    }

    /**
     * 构造函数
     *
     * @param lon     经度
     * @param lat     纬度
     * @param ts      时间戳
     * @param speed   速度
     * @param heading 方向角
     */
    public TrackPoint(long ts, double lon, double lat, double speed, int heading) {
        super(lon, lat);
        this.ts = ts;
        this.speed = speed;
        this.heading = heading;
    }

    /**
     * 构造函数，创建一个轨迹点对象
     *
     * @param lon 经度坐标值
     * @param lat 纬度坐标值
     */
    public TrackPoint(double lon, double lat) {
        super(lon, lat);
    }

    public long getTs() {
        return ts;
    }

    public void setTs(long ts) {
        this.ts = ts;
    }

    public double getSpeed() {
        return speed;
    }

    public void setSpeed(double speed) {
        this.speed = speed;
    }

    public boolean getValid() {
        return valid;
    }

    public void setValid(boolean valid) {
        this.valid = valid;
    }


    public int getHeading() {
        return heading;
    }

    public void setHeading(int heading) {
        this.heading = heading;
    }
}
