package sunyu.util.pojo;

import java.time.LocalDateTime;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;

import org.locationtech.jts.geom.Geometry;

import cn.hutool.core.collection.CollUtil;
import cn.hutool.core.date.LocalDateTimeUtil;
import cn.hutool.core.util.StrUtil;

/**
 * 轮廓区块详情：单个 Polygon 的元数据。
 * 
 * @author SunYu
 */
public class OutlinePart {
    /**
     * 区块的几何图形（Polygon）
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
     * 区块的WKT表示
     */
    private String wkt;
    /**
     * 轮廓内的轨迹点集合（按时间升序）
     */
    private List<TrackPoint> trackPoints;
    /**
     * 轮廓内的轨迹点字符串表示（按时间升序）
     */
    private String trackStr;
    /**
     * 使用作业宽幅（米）
     */
    private double totalWidthM;

    public OutlinePart(Geometry outline, LocalDateTime startTime, LocalDateTime endTime, double mu, String wkt) {
        this.outline = outline;
        this.startTime = startTime;
        this.endTime = endTime;
        this.mu = mu;
        this.wkt = wkt;
    }

    public OutlinePart(Geometry outline, LocalDateTime startTime, LocalDateTime endTime, double mu, String wkt,
            List<TrackPoint> trackPoints) {
        this(outline, startTime, endTime, mu, wkt);
        if (trackPoints != null) {
            List<TrackPoint> copy = new ArrayList<>(trackPoints);
            copy.sort(Comparator.comparing(
                    TrackPoint::getTime,
                    Comparator.nullsLast(Comparator.naturalOrder())));
            this.trackPoints = copy;
        } else {
            this.trackPoints = null;
        }
    }

    public Geometry getOutline() {
        return outline;
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

    public List<TrackPoint> getTrackPoints() {
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

    public double getTotalWidthM() {
        return totalWidthM;
    }

    public OutlinePart setTotalWidthM(double totalWidthM) {
        this.totalWidthM = totalWidthM;
        return this;
    }
}
