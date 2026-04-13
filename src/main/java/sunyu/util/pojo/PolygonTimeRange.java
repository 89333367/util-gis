package sunyu.util.pojo;

import java.time.LocalDateTime;

public class PolygonTimeRange {
    private Integer polygonIndex;
    private LocalDateTime start;
    private LocalDateTime end;

    public PolygonTimeRange(Integer polygonIndex, LocalDateTime start, LocalDateTime end) {
        this.polygonIndex = polygonIndex;
        this.start = start;
        this.end = end;
    }

    public Integer getPolygonIndex() {
        return polygonIndex;
    }

    public LocalDateTime getStart() {
        return start;
    }

    public LocalDateTime getEnd() {
        return end;
    }
}
