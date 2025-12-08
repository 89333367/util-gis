package sunyu.util.pojo;

import java.time.LocalDateTime;
import java.util.List;

public class Part {
    /**
     * 区块的起始时间
     */
    private LocalDateTime startTime;
    /**
     * 区块的结束时间
     */
    private LocalDateTime endTime;
    /**
     * 区块面积（亩）
     */
    private double mu;
    /**
     * 区块的WKT表示（WGS84坐标系）
     */
    private String wkt;
    /**
     * 区块内的轨迹点集合（WGS84坐标系）
     */
    private List<Wgs84Point> trackPoints;
    /**
     * 区块内的轨迹点字符串表示形式（经度,纬度,定位时间#经度,纬度,定位时间#）
     * <p>
     * 定位时间格式yyyyMMddHHmmss
     */
    private String trackStr;

    public LocalDateTime getStartTime() {
        return startTime;
    }

    public void setStartTime(LocalDateTime startTime) {
        this.startTime = startTime;
    }

    public LocalDateTime getEndTime() {
        return endTime;
    }

    public void setEndTime(LocalDateTime endTime) {
        this.endTime = endTime;
    }

    public double getMu() {
        return mu;
    }

    public void setMu(double mu) {
        this.mu = mu;
    }

    public String getWkt() {
        return wkt;
    }

    public void setWkt(String wkt) {
        this.wkt = wkt;
    }

    public List<Wgs84Point> getTrackPoints() {
        return trackPoints;
    }

    public void setTrackPoints(List<Wgs84Point> trackPoints) {
        this.trackPoints = trackPoints;
    }

    public String getTrackStr() {
        return trackStr;
    }

    public void setTrackStr(String trackStr) {
        this.trackStr = trackStr;
    }

}
