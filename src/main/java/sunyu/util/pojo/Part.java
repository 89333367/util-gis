package sunyu.util.pojo;

import cn.hutool.core.collection.CollUtil;
import cn.hutool.core.date.LocalDateTimeUtil;
import cn.hutool.core.util.StrUtil;

import java.time.LocalDateTime;
import java.util.Comparator;
import java.util.List;

/**
 * 地块信息
 */
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
        trackPoints.sort(Comparator.comparing(Wgs84Point::getGpsTime));
        this.trackPoints = trackPoints;
    }

    public String getTrackStr() {
        if (CollUtil.isNotEmpty(trackPoints)) {
            StringBuilder sb = new StringBuilder();
            for (Wgs84Point trackPoints : trackPoints) {
                sb.append(StrUtil.format("{},{},{}#", trackPoints.getLongitude(), trackPoints.getLatitude(), LocalDateTimeUtil.format(trackPoints.getGpsTime(), "yyyyMMddHHmmss")));
            }
            return sb.toString();
        }
        return "";
    }

    public void setTrackStr(String trackStr) {
        this.trackStr = trackStr;
    }

}
