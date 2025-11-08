package sunyu.util.pojo;

import java.time.LocalDateTime;
import java.util.List;

import org.locationtech.jts.geom.Geometry;

/**
 * 分割道路的结果，包含整体轮廓以及每个区块的详情
 * 
 * @author SunYu
 */
public class SplitRoadResult {
    /**
     * 整体轮廓的几何图形（高斯投影）
     */
    private Geometry outline;
    /**
     * 整体轮廓的起始时间
     */
    private LocalDateTime startTime;
    /**
     * 整体轮廓的结束时间
     */
    private LocalDateTime endTime;
    /**
     * 每个区块的详情
     */
    private List<OutlinePart> parts;
    /**
     * 轮廓的WKT表示（WGS84坐标系）
     */
    private String wkt;
    /**
     * 使用作业总宽幅（米）
     */
    private double totalWidthM;

    /**
     * 整体轮廓面积（亩）
     */
    private double mu;

    public SplitRoadResult() {
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

    public List<OutlinePart> getParts() {
        return parts;
    }

    public void setParts(List<OutlinePart> parts) {
        this.parts = parts;
    }

    public String getWkt() {
        return wkt;
    }

    public void setWkt(String wkt) {
        this.wkt = wkt;
    }

    public double getTotalWidthM() {
        return totalWidthM;
    }

    public void setTotalWidthM(double totalWidthM) {
        this.totalWidthM = totalWidthM;
    }

    public double getMu() {
        return mu;
    }

    public void setMu(double mu) {
        this.mu = mu;
    }

}