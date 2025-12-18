package sunyu.util.pojo;

import org.locationtech.jts.geom.Geometry;

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
    private List<SplitPart> splitParts = new ArrayList<>();

    /**
     * 地块的高斯投影几何图形（高斯投影坐标系）
     */
    private Geometry gaussGeometry;

    /**
     * 最小有效时间间隔（秒）
     */
    private int minEffectiveInterval;

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

    public List<SplitPart> getParts() {
        splitParts.sort(Comparator.comparing(SplitPart::getStartTime));
        return splitParts;
    }

    public void setParts(List<SplitPart> splitParts) {
        this.splitParts = splitParts;
    }

    public Geometry getGaussGeometry() {
        return gaussGeometry;
    }

    public void setGaussGeometry(Geometry gaussGeometry) {
        this.gaussGeometry = gaussGeometry;
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

    public List<SplitPart> getSplitParts() {
        return splitParts;
    }

    public void setSplitParts(List<SplitPart> splitParts) {
        this.splitParts = splitParts;
    }

    public int getMinEffectiveInterval() {
        return minEffectiveInterval;
    }

    public void setMinEffectiveInterval(int minEffectiveInterval) {
        this.minEffectiveInterval = minEffectiveInterval;
    }
}
