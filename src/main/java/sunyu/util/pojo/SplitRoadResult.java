package sunyu.util.pojo;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;

import org.locationtech.jts.geom.Geometry;

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

    public double getTotalWidthM() {
        return totalWidthM;
    }

    public SplitRoadResult setTotalWidthM(double totalWidthM) {
        this.totalWidthM = totalWidthM;
        return this;
    }
}