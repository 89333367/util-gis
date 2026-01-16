package sunyu.util.pojo;

/**
 * 道路切割参数
 *
 * @author SunYu
 */
public class SplitRoadParams {
    /**
     * DBSCAN半径（米）
     * <p>
     * 用于DBSCAN聚类算法的半径阈值。
     * 当两个点之间的距离小于此值时，被认为是邻居点，属于同一个聚类。
     * </p>
     */
    private Double dbScanEpsilon;
    /**
     * DBSCAN最小点数
     * <p>
     * 用于DBSCAN聚类算法的最小点阈值。
     * 当一个区域内的点数量小于此值时，被认为是噪声点或异常值，不会被分配到任何聚类中。
     * </p>
     */
    private Integer dbScanMinPoints;
    /**
     * 道路宽度(米)
     * <p>
     * 用于定义切割后的道路宽度。
     * 切割时，会根据此宽度将道路切割掉。
     * </p>
     */
    private Double roadWidth;
    /**
     * 最小返回多边形面积限制(亩)
     * <p>
     * 小于此面积的多边形会被过滤掉
     * </p>
     */
    private Double minReturnMu;
    /**
     * 是否检查作业状态
     * <p>
     * 如果为true，则会根据作业状态来判断是否切割。
     * 如果为false，则会忽略作业状态，直接切割。
     * </p>
     */
    private boolean checkWorkingStatus = true;

    public SplitRoadParams() {
    }

    /**
     * 构造函数
     *
     * @param roadWidth 道路宽度(米)
     */
    public SplitRoadParams(Double roadWidth) {
        this.roadWidth = roadWidth;
    }

    /**
     * 构造函数
     *
     * @param dbScanEpsilon   DBSCAN半径（米）
     * @param dbScanMinPoints DBSCAN最小点数
     */
    public SplitRoadParams(Double dbScanEpsilon, Integer dbScanMinPoints) {
        this.dbScanEpsilon = dbScanEpsilon;
        this.dbScanMinPoints = dbScanMinPoints;
    }

    /**
     * 构造函数
     *
     * @param dbScanEpsilon   DBSCAN半径（米）
     * @param dbScanMinPoints DBSCAN最小点数
     * @param roadWidth       道路宽度(米)
     */
    public SplitRoadParams(Double dbScanEpsilon, Integer dbScanMinPoints, Double roadWidth) {
        this.dbScanEpsilon = dbScanEpsilon;
        this.dbScanMinPoints = dbScanMinPoints;
        this.roadWidth = roadWidth;
    }

    /**
     * 构造函数
     *
     * @param dbScanEpsilon   DBSCAN半径（米）
     * @param dbScanMinPoints DBSCAN最小点数
     * @param roadWidth       道路宽度(米)
     * @param minReturnMu     最小返回多边形面积限制(亩)
     */
    public SplitRoadParams(Double dbScanEpsilon, Integer dbScanMinPoints, Double roadWidth, Double minReturnMu) {
        this.dbScanEpsilon = dbScanEpsilon;
        this.dbScanMinPoints = dbScanMinPoints;
        this.roadWidth = roadWidth;
        this.minReturnMu = minReturnMu;
    }

    /**
     * 构造函数
     *
     * @param dbScanEpsilon      DBSCAN半径（米）
     * @param dbScanMinPoints    DBSCAN最小点数
     * @param roadWidth          道路宽度(米)
     * @param minReturnMu        最小返回多边形面积限制(亩)
     * @param checkWorkingStatus 是否检查作业状态
     */
    public SplitRoadParams(Double dbScanEpsilon, Integer dbScanMinPoints, Double roadWidth, Double minReturnMu, boolean checkWorkingStatus) {
        this.dbScanEpsilon = dbScanEpsilon;
        this.dbScanMinPoints = dbScanMinPoints;
        this.roadWidth = roadWidth;
        this.minReturnMu = minReturnMu;
        this.checkWorkingStatus = checkWorkingStatus;
    }

    public Double getDbScanEpsilon() {
        return dbScanEpsilon;
    }

    public void setDbScanEpsilon(Double dbScanEpsilon) {
        this.dbScanEpsilon = dbScanEpsilon;
    }

    public Integer getDbScanMinPoints() {
        return dbScanMinPoints;
    }

    public void setDbScanMinPoints(Integer dbScanMinPoints) {
        this.dbScanMinPoints = dbScanMinPoints;
    }

    public Double getRoadWidth() {
        return roadWidth;
    }

    public void setRoadWidth(Double roadWidth) {
        this.roadWidth = roadWidth;
    }

    public Double getMinReturnMu() {
        return minReturnMu;
    }

    public void setMinReturnMu(Double minReturnMu) {
        this.minReturnMu = minReturnMu;
    }

    public boolean getCheckWorkingStatus() {
        return checkWorkingStatus;
    }

    public void setCheckWorkingStatus(boolean checkWorkingStatus) {
        this.checkWorkingStatus = checkWorkingStatus;
    }
}