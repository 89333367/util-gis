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