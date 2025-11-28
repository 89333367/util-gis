package sunyu.util.pojo;

import org.apache.commons.math3.ml.clustering.Clusterable;

import java.time.LocalDateTime;

/**
 * 轨迹点数据类，继承自CoordinatePoint，实现Clusterable接口以支持Apache Commons Math聚类
 *
 * @author SunYu
 */
public class TrackPoint extends CoordinatePoint implements Clusterable {
    /**
     * 定位时间
     */
    private LocalDateTime time;

    /**
     * 两点间的距离(米)
     */
    private double distance;

    /**
     * 两点间的速度(米/秒)
     */
    private double speed;

    /**
     * 两点间的方向角(度)
     */
    private double direction;

    public TrackPoint() {
    }

    /**
     * 完整构造函数：包含所有属性
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

    /**
     * 两点间的距离(米)
     */
    public double getDistance() {
        return distance;
    }

    /**
     * 两点间的距离(米)
     */
    public void setDistance(double distance) {
        this.distance = distance;
    }

    /**
     * 两点间的速度(米/秒)
     */
    public double getSpeed() {
        return speed;
    }

    /**
     * 两点间的速度(米/秒)
     */
    public void setSpeed(double speed) {
        this.speed = speed;
    }

    /**
     * 两点间的方向角(度)
     */
    public double getDirection() {
        return direction;
    }

    /**
     * 两点间的方向角(度)
     */
    public void setDirection(double direction) {
        this.direction = direction;
    }

    /**
     * 实现Clusterable接口，返回用于聚类的点坐标
     *
     * @return 包含经度和纬度的双精度数组
     */
    @Override
    public double[] getPoint() {
        return new double[]{getLon(), getLat()};
    }

}