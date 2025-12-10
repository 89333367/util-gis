package sunyu.util.pojo;

import java.time.LocalDateTime;

/**
 * 高斯投影点
 */
public class GaussPoint extends Wgs84Point {
    /**
     * 高斯投影x坐标(单位：米)，内部使用，外部不要直接访问
     */
    private double gaussX;
    /**
     * 高斯投影y坐标(单位：米)，内部使用，外部不要直接访问
     */
    private double gaussY;

    public GaussPoint() {
    }

    /**
     * 构造方法
     *
     * @param gpsTime   GPS定位时间
     * @param longitude 经度(单位：度)(东经为正值，西经为负值)(范围：-180.0到180.0)
     * @param latitude  纬度(单位：度)(北纬为正值，南纬为负值)(范围：-90.0到90.0)
     * @param gaussX    高斯投影x坐标(单位：米)
     * @param gaussY    高斯投影y坐标(单位：米)
     */
    public GaussPoint(LocalDateTime gpsTime, double longitude, double latitude, double gaussX, double gaussY) {
        super(gpsTime, longitude, latitude);
        this.gaussX = gaussX;
        this.gaussY = gaussY;
    }

    public double getGaussX() {
        return gaussX;
    }

    public void setGaussX(double gaussX) {
        this.gaussX = gaussX;
    }

    public double getGaussY() {
        return gaussY;
    }

    public void setGaussY(double gaussY) {
        this.gaussY = gaussY;
    }
}