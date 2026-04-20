package sunyu.util.pojo;

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

    /**
     * 点位属于哪个多边形
     */
    private Integer polygonIndex;

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

    public Integer getPolygonIndex() {
        return polygonIndex;
    }

    public void setPolygonIndex(Integer polygonIndex) {
        this.polygonIndex = polygonIndex;
    }
}