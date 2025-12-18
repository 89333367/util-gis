package sunyu.util.pojo;

import org.locationtech.jts.geom.Geometry;

import java.time.LocalDateTime;

/**
 * 地块信息
 */
public class SplitPart {
    /**
     * 地块的起始时间
     */
    private LocalDateTime startTime;
    /**
     * 地块的结束时间
     */
    private LocalDateTime endTime;
    /**
     * 地块面积（亩）
     */
    private double mu;
    /**
     * 地块的WKT表示（WGS84坐标系）
     */
    private String wkt;
    /**
     * 地块的高斯投影几何图形（高斯投影坐标系）
     */
    private Geometry gaussGeometry;

    /**
     * 最小有效时间间隔（秒）
     */
    private int minEffectiveInterval;

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

    public Geometry getGaussGeometry() {
        return gaussGeometry;
    }

    public void setGaussGeometry(Geometry gaussGeometry) {
        this.gaussGeometry = gaussGeometry;
    }

    public int getMinEffectiveInterval() {
        return minEffectiveInterval;
    }

    public void setMinEffectiveInterval(int minEffectiveInterval) {
        this.minEffectiveInterval = minEffectiveInterval;
    }
}
