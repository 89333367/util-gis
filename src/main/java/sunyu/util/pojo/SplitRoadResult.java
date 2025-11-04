package sunyu.util.pojo;

import java.time.LocalDateTime;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;

import org.locationtech.jts.geom.Geometry;

import cn.hutool.core.collection.CollUtil;

/**
 * 分割道路的结果，包含整体轮廓以及每个区块的详情（起止时间、亩数）
 * 
 * @author SunYu
 */
public class SplitRoadResult {
    /**
     * 整体轮廓的几何图形
     */
    private final Geometry outline;
    /**
     * 区块的起始时间
     */
    private LocalDateTime startTime;
    /**
     * 区块的结束时间
     */
    private LocalDateTime endTime;
    /**
     * 每个区块的详情
     */
    private final List<OutlinePart> parts;
    /**
     * 轮廓的WKT表示
     */
    private final String wkt;
    /**
     * 使用作业宽幅（米）
     */
    private double totalWidthM;

    /**
     * 区块面积（亩）
     */
    private double mu;

    public SplitRoadResult(Geometry outline, List<OutlinePart> parts) {
        this(outline, parts, null);
    }

    public SplitRoadResult(Geometry outline, List<OutlinePart> parts, String wkt) {
        this.outline = outline;
        this.parts = parts;
        this.wkt = wkt;
    }

    public Geometry getOutline() {
        return outline;
    }

    public List<OutlinePart> getParts() {
        if (parts == null || parts.isEmpty()) {
            return parts;
        }
        List<OutlinePart> sorted = new ArrayList<>(parts);
        sorted.sort(
                Comparator.comparing(
                        OutlinePart::getStartTime,
                        Comparator.nullsLast(Comparator.naturalOrder())));
        return sorted;
    }

    public String getWkt() {
        return wkt;
    }

    public double getMu() {
        double mu = 0.0;
        if (CollUtil.isNotEmpty(parts)) {
            for (OutlinePart part : parts) {
                mu += part.getMu();
            }
        }
        // 返回保留4位小数的精度
        return Math.round(mu * 10000.0) / 10000.0;
    }

    public SplitRoadResult setTotalWidthM(double totalWidthM) {
        this.totalWidthM = totalWidthM;
        return this;
    }

    public double getTotalWidthM() {
        return totalWidthM;
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

    /**
     * 根据parts列表计算并设置整体的起止时间
     */
    public void calculateTimeRange() {
        if (parts == null || parts.isEmpty()) {
            return;
        }

        LocalDateTime earliestStart = null;
        LocalDateTime latestEnd = null;

        for (OutlinePart part : parts) {
            if (part.getStartTime() != null) {
                if (earliestStart == null || part.getStartTime().isBefore(earliestStart)) {
                    earliestStart = part.getStartTime();
                }
            }

            if (part.getEndTime() != null) {
                if (latestEnd == null || part.getEndTime().isAfter(latestEnd)) {
                    latestEnd = part.getEndTime();
                }
            }
        }

        this.startTime = earliestStart;
        this.endTime = latestEnd;
    }

}