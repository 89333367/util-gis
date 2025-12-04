package sunyu.util.pojo;

import java.time.LocalDateTime;

public class Wgs84Point {
    /**
     * GPS定位时间
     */
    private LocalDateTime gpsTime;
    /**
     * 经度(单位：度)(东经为正值，西经为负值)(范围：-180.0到180.0)
     */
    private double longitude;
    /**
     * 纬度(单位：度)(北纬为正值，南纬为负值)(范围：-90.0到90.0)
     */
    private double latitude;
    /**
     * GPS定位状态(0：未知；1：已定位；2：未定位)
     */
    private int gpsStatus;
    /**
     * 作业状态(0：未知；1：作业中；2：未作业)
     */
    private int jobStatus;

    public Wgs84Point() {
    }

    /**
     * 构造方法
     *
     * @param gpsTime   GPS定位时间
     * @param longitude 经度(单位：度)(东经为正值，西经为负值)(范围：-180.0到180.0)
     * @param latitude  纬度(单位：度)(北纬为正值，南纬为负值)(范围：-90.0到90.0)
     */
    public Wgs84Point(LocalDateTime gpsTime, double longitude, double latitude) {
        this.gpsTime = gpsTime;
        this.longitude = longitude;
        this.latitude = latitude;
    }

    /**
     * 构造方法
     *
     * @param gpsTime   GPS定位时间
     * @param longitude 经度(单位：度)(东经为正值，西经为负值)(范围：-180.0到180.0)
     * @param latitude  纬度(单位：度)(北纬为正值，南纬为负值)(范围：-90.0到90.0)
     * @param gpsStatus GPS定位状态(0：未知；1：已定位；2：未定位)
     */
    public Wgs84Point(LocalDateTime gpsTime, double longitude, double latitude, int gpsStatus) {
        this.gpsTime = gpsTime;
        this.longitude = longitude;
        this.latitude = latitude;
        this.gpsStatus = gpsStatus;
    }

    /**
     * 构造方法
     *
     * @param gpsTime   GPS定位时间
     * @param longitude 经度(单位：度)(东经为正值，西经为负值)(范围：-180.0到180.0)
     * @param latitude  纬度(单位：度)(北纬为正值，南纬为负值)(范围：-90.0到90.0)
     * @param gpsStatus GPS定位状态(0：未知；1：已定位；2：未定位)
     * @param jobStatus 作业状态(0：未知；1：作业中；2：未作业)
     */
    public Wgs84Point(LocalDateTime gpsTime, double longitude, double latitude, int gpsStatus, int jobStatus) {
        this.gpsTime = gpsTime;
        this.longitude = longitude;
        this.latitude = latitude;
        this.gpsStatus = gpsStatus;
        this.jobStatus = jobStatus;
    }

    public LocalDateTime getGpsTime() {
        return gpsTime;
    }

    public void setGpsTime(LocalDateTime gpsTime) {
        this.gpsTime = gpsTime;
    }

    public double getLongitude() {
        return longitude;
    }

    public void setLongitude(double longitude) {
        this.longitude = longitude;
    }

    public double getLatitude() {
        return latitude;
    }

    public void setLatitude(double latitude) {
        this.latitude = latitude;
    }

    public int getGpsStatus() {
        return gpsStatus;
    }

    public void setGpsStatus(int gpsStatus) {
        this.gpsStatus = gpsStatus;
    }

    public int getJobStatus() {
        return jobStatus;
    }

    public void setJobStatus(int jobStatus) {
        this.jobStatus = jobStatus;
    }

}