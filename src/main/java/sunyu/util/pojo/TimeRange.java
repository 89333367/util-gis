package sunyu.util.pojo;

import java.time.LocalDateTime;
import java.util.List;

public class TimeRange {
    private LocalDateTime start;
    private LocalDateTime end;
    private List<GaussPoint> gaussPoints;

    public TimeRange(LocalDateTime start, LocalDateTime end, List<GaussPoint> gaussPoints) {
        this.start = start;
        this.end = end;
        this.gaussPoints = gaussPoints;
    }

    public TimeRange(LocalDateTime start, LocalDateTime end) {
        this.start = start;
        this.end = end;
    }

    public TimeRange() {
    }

    public LocalDateTime getStart() {
        return start;
    }

    public void setStart(LocalDateTime start) {
        this.start = start;
    }

    public LocalDateTime getEnd() {
        return end;
    }

    public void setEnd(LocalDateTime end) {
        this.end = end;
    }

    public List<GaussPoint> getGaussPoints() {
        return gaussPoints;
    }

    public void setGaussPoints(List<GaussPoint> gaussPoints) {
        this.gaussPoints = gaussPoints;
    }
}
