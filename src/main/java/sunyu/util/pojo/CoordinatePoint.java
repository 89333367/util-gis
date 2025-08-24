package sunyu.util.pojo;

/**
 * 坐标点基类，用于存储经纬度坐标对
 *
 * @author SunYu
 */
public class CoordinatePoint {
    /**
     * 经度坐标（东经为正，西经为负）
     */
    private double lon;

    /**
     * 纬度坐标（北纬为正，南纬为负）
     */
    private double lat;

    /**
     * 构造函数
     *
     * @param lon 经度
     * @param lat 纬度
     */
    public CoordinatePoint(double lon, double lat) {
        this.lon = lon;
        this.lat = lat;
    }

    public double getLon() {
        return lon;
    }

    public void setLon(double lon) {
        this.lon = lon;
    }

    public double getLat() {
        return lat;
    }

    public void setLat(double lat) {
        this.lat = lat;
    }
}
