package sunyu.util.pojo;

import java.time.LocalDateTime;

import org.locationtech.jts.geom.Geometry;

import cn.hutool.core.collection.CollUtil;
import cn.hutool.core.date.LocalDateTimeUtil;
import cn.hutool.core.util.StrUtil;

/**
 * 轮廓区块详情：单个 Polygon 的元数据。
 * @author SunYu
 */
public class OutlinePart {
    /**
     * 区块的几何图形（Polygon）
     */
    private Geometry polygon;
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
     * 区块的WKT表示
     */
    private String wkt;
    /**
     * 轮廓内的轨迹点集合（按时间升序）
     */
    private java.util.List<TrackPoint> trackPoints;

    private String trackStr;

    public OutlinePart(Geometry polygon, LocalDateTime startTime, LocalDateTime endTime, double mu, String wkt) {
        this.polygon = polygon;
        this.startTime = startTime;
        this.endTime = endTime;
        this.mu = mu;
        this.wkt = wkt;
    }

    public OutlinePart(Geometry polygon, LocalDateTime startTime, LocalDateTime endTime, double mu, String wkt,
            java.util.List<TrackPoint> trackPoints) {
        this(polygon, startTime, endTime, mu, wkt);
        if (trackPoints != null) {
            java.util.List<TrackPoint> copy = new java.util.ArrayList<>(trackPoints);
            copy.sort(java.util.Comparator.comparing(
                    TrackPoint::getTime,
                    java.util.Comparator.nullsLast(java.util.Comparator.naturalOrder())));
            this.trackPoints = copy;
        } else {
            this.trackPoints = null;
        }
    }

    public Geometry getPolygon() {
        return polygon;
    }

    public LocalDateTime getStartTime() {
        return startTime;
    }

    public LocalDateTime getEndTime() {
        return endTime;
    }

    public double getMu() {
        return mu;
    }

    public String getWkt() {
        return wkt;
    }

    public java.util.List<TrackPoint> getTrackPoints() {
        return trackPoints;
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

}
