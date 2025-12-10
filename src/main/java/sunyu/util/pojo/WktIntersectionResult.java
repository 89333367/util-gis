package sunyu.util.pojo;

/**
 * 相交地块结果类
 */
public class WktIntersectionResult {
    /**
     * 相交地块WKT(WGS84坐标系)
     */
    private String wkt;
    /**
     * 相交的面积（亩）
     */
    private double mu;

    public WktIntersectionResult() {
    }

    public WktIntersectionResult(String wkt, double mu) {
        this.wkt = wkt;
        this.mu = mu;
    }

    public String getWkt() {
        return wkt;
    }

    public void setWkt(String wkt) {
        this.wkt = wkt;
    }

    public double getMu() {
        return mu;
    }

    public void setMu(double mu) {
        this.mu = mu;
    }

}
