package sunyu.util.pojo;

import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.Point;

import java.time.LocalDateTime;

/**
 * 地块信息
 */
public class FarmPlot {
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
     * 地块的WGS84坐标系几何图形
     */
    private Geometry wgs84Geometry;
    /**
     * 最小有效时间间隔（秒）
     */
    private int minEffectiveInterval;
    /**
     * 聚类点的数量
     */
    private int clusterPointCount;
    /**
     * 使用作业总宽幅（米）
     */
    private double workingWidth;
    /**
     * 中心点（WGS84坐标系）
     */
    private Wgs84Point centerWgs84Point;

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


    public double getWorkingWidth() {
        return workingWidth;
    }

    public void setWorkingWidth(double workingWidth) {
        this.workingWidth = workingWidth;
    }

    public int getClusterPointCount() {
        return clusterPointCount;
    }

    public void setClusterPointCount(int clusterPointCount) {
        this.clusterPointCount = clusterPointCount;
    }

    public Geometry getWgs84Geometry() {
        return wgs84Geometry;
    }

    public void setWgs84Geometry(Geometry wgs84Geometry) {
        this.wgs84Geometry = wgs84Geometry;
    }

    public Wgs84Point getCenterWgs84Point() {
        if (wgs84Geometry != null) {
            Wgs84Point wgs84Point = new Wgs84Point();
            Point centerPoint = wgs84Geometry.getInteriorPoint();
            wgs84Point.setLongitude(centerPoint.getX());
            wgs84Point.setLatitude(centerPoint.getY());
            return wgs84Point;
        }
        return centerWgs84Point;
    }

    public void setCenterWgs84Point(Wgs84Point centerWgs84Point) {
        this.centerWgs84Point = centerWgs84Point;
    }
}