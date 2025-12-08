package sunyu.util.pojo;

/**
 * 相交轮廓结果类
 */
public class WktIntersectionResult {
    /**
     * 相交轮廓WKT(WGS84坐标系)
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
