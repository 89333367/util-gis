package sunyu.util.pojo;

import java.time.LocalDateTime;
import java.util.List;

import org.locationtech.jts.geom.Geometry;

import cn.hutool.core.collection.CollUtil;
import cn.hutool.core.date.LocalDateTimeUtil;
import cn.hutool.core.util.StrUtil;

/**
 * 单个轮廓区块详情
 * 
 * @author SunYu
 */
public class OutlinePart {
    /**
     * 区块的几何图形（Polygon，平面投影，可能是高斯投影或者是UTM投影）
     */
    private Geometry outline;
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
    private List<TrackPoint> trackPoints;
    /**
     * 区块内的轨迹点字符串表示形式（经度,纬度,定位时间#经度,纬度,定位时间#）
     * 
     * 定位时间格式yyyyMMddHHmmss
     */
    private String trackStr;
    /**
     * 使用作业总宽幅（米）
     */
    private double totalWidthM;

    public OutlinePart() {
    }

    public Geometry getOutline() {
        return outline;
    }

    public void setOutline(Geometry outline) {
        this.outline = outline;
    }

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

    public List<TrackPoint> getTrackPoints() {
        return trackPoints;
    }

    public void setTrackPoints(List<TrackPoint> trackPoints) {
        if (CollUtil.isNotEmpty(trackPoints)) {
            trackPoints.sort((a, b) -> a.getTime().compareTo(b.getTime()));
        }
        this.trackPoints = trackPoints;
    }

    public String getTrackStr() {
        if (CollUtil.isNotEmpty(trackPoints)) {
            StringBuilder sb = new StringBuilder();
            for (TrackPoint trackPoints : trackPoints) {
                sb.append(StrUtil.format("{},{},{}#", trackPoints.getLon(), trackPoints.getLat(),
                        LocalDateTimeUtil.format(trackPoints.getTime(), "yyyyMMddHHmmss")));
            }
            trackStr = sb.toString();
        }
        return trackStr;
    }

    public void setTrackStr(String trackStr) {
        this.trackStr = trackStr;
    }

    public double getTotalWidthM() {
        return totalWidthM;
    }

    public void setTotalWidthM(double totalWidthM) {
        this.totalWidthM = totalWidthM;
    }

}
