package sunyu.util.pojo;

import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.Point;

import java.time.LocalDateTime;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;

/**
 * 道路拆分结果
 */
public class SplitResult {
    /**
     * 地块的起始时间
     */
    private LocalDateTime startTime;
    /**
     * 地块的结束时间
     */
    private LocalDateTime endTime;
    /**
     * 轮廓的WKT表示（WGS84坐标系）
     */
    private String wkt;
    /**
     * 使用作业总宽幅（米）
     */
    private double workingWidth;
    /**
     * 整体轮廓面积（亩）
     */
    private double mu;
    /**
     * 拆分后的地块列表
     */
    private List<FarmPlot> farmPlots = new ArrayList<>();
    /**
     * 地块的WGS84坐标系几何图形
     */
    private Geometry wgs84Geometry;
    /**
     * 聚类点的数量
     */
    private int clusterPointCount;
    /**
     * 中心点（WGS84坐标系）
     */
    private Wgs84Point centerWgs84Point;

    public String getWkt() {
        return wkt;
    }

    public void setWkt(String wkt) {
        this.wkt = wkt;
    }

    public double getWorkingWidth() {
        return workingWidth;
    }

    public void setWorkingWidth(double workingWidth) {
        this.workingWidth = workingWidth;
    }

    public double getMu() {
        return mu;
    }

    public void setMu(double mu) {
        this.mu = mu;
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

    public void setSplitParts(List<FarmPlot> farmPlots) {
        this.farmPlots = farmPlots;
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

    public int getClusterPointCount() {
        // 返回总聚类点数
        return farmPlots.stream().mapToInt(FarmPlot::getClusterPointCount).sum();
    }

    public void setClusterPointCount(int clusterPointCount) {
        this.clusterPointCount = clusterPointCount;
    }

    public List<FarmPlot> getFarmPlots() {
        // 返回保持作业时间升序
        farmPlots.sort(Comparator.comparing(FarmPlot::getStartTime));
        return farmPlots;
    }

    public void setFarmPlots(List<FarmPlot> farmPlots) {
        this.farmPlots = farmPlots;
    }

}