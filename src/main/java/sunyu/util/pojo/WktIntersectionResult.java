package sunyu.util.pojo;

public class WktIntersectionResult {
    private String wkt;
    private double mu;

    public WktIntersectionResult() {}

    public WktIntersectionResult(String wkt, double mu) {
        this.wkt = wkt;
        this.mu = mu;
    }

    public String getWkt() {
        return wkt;
    }

    public WktIntersectionResult setWkt(String wkt) {
        this.wkt = wkt;
        return this;
    }

    public double getMu() {
        return mu;
    }

    public WktIntersectionResult setMu(double mu) {
        this.mu = mu;
        return this;
    }

    @Override
    public String toString() {
        return "WktIntersectionResult{" +
                "wkt='" + wkt + '\'' +
                ", mu=" + mu +
                '}';
    }
}
