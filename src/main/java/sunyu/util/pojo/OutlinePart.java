package sunyu.util.pojo;

import java.time.LocalDateTime;

import org.locationtech.jts.geom.Geometry;

/**
 * 轮廓区块详情：单个 Polygon 的元数据。
 * @author SunYu
 */
public class OutlinePart {
    /**
     * 区块的几何图形（Polygon）
     */
    private final Geometry polygon;
    /**
     * 区块的起始时间
     */
    private final LocalDateTime startTime;
    /**
     * 区块的结束时间
     */
    private final LocalDateTime endTime;
    /**
     * 区块面积（亩）
     */
    private final double mu;
    /**
     * 区块的WKT表示
     */
    private final String wkt;

    public OutlinePart(Geometry polygon, LocalDateTime startTime, LocalDateTime endTime, double mu, String wkt) {
        this.polygon = polygon;
        this.startTime = startTime;
        this.endTime = endTime;
        this.mu = mu;
        this.wkt = wkt;
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
}
