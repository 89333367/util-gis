package sunyu.util.pojo;

import java.util.List;

import org.locationtech.jts.geom.Geometry;

/**
 * 分割道路的结果，包含整体轮廓以及每个区块的详情（起止时间、亩数）
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
        java.util.List<OutlinePart> sorted = new java.util.ArrayList<>(parts);
        sorted.sort(
                java.util.Comparator.comparing(
                        OutlinePart::getStartTime,
                        java.util.Comparator.nullsLast(java.util.Comparator.naturalOrder())));
        return sorted;
    }

    public String getWkt() {
        return wkt;
    }
}
