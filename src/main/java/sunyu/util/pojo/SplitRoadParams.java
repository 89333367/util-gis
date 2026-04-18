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
     * 正缓冲
     */
    private Double positiveBuffer;
    /**
     * 负缓冲
     */
    private Double negativeBuffer;
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
    private Boolean checkWorkingStatus = false;

    /**
     * 算法索引
     * 0：代表如果地块之间有时间交叉，那么将时间交叉的地块进行合并
     * 1：代表如果地块之间有时间交叉，那么将地块内再次进行时间窗口拆分
     */
    private int algorithmIndex = 1;

    public SplitRoadParams() {
    }

    public Double getDbScanEpsilon() {
        return dbScanEpsilon;
    }

    /**
     * 聚类半径（米）
     * <p>
     * 用于空间密集聚类算法的半径阈值。
     * 当两个点之间的距离小于此值时，被认为是邻居点，属于同一个聚类。
     * </p>
     */
    public SplitRoadParams setDbScanEpsilon(Double dbScanEpsilon) {
        this.dbScanEpsilon = dbScanEpsilon;
        return this;
    }

    public Integer getDbScanMinPoints() {
        return dbScanMinPoints;
    }

    /**
     * 聚类最小点数设置
     * <p>
     * 用于空间密集聚类算法的最小点阈值。
     * 当一个区域内的点数量小于此值时，被认为是噪声点或异常值，不会被分配到任何聚类中。
     * </p>
     */
    public SplitRoadParams setDbScanMinPoints(Integer dbScanMinPoints) {
        this.dbScanMinPoints = dbScanMinPoints;
        return this;
    }

    public Double getMinReturnMu() {
        return minReturnMu;
    }

    /**
     * 最小返回多边形面积限制(亩)
     * <p>
     * 小于此面积的多边形会被过滤掉
     * </p>
     */
    public SplitRoadParams setMinReturnMu(Double minReturnMu) {
        this.minReturnMu = minReturnMu;
        return this;
    }

    public Boolean getCheckWorkingStatus() {
        return checkWorkingStatus;
    }

    /**
     * 是否检查作业状态
     * <p>
     * 如果为true，则抛弃未作业状态的点
     * 如果为false，则聚类计算时会计算未作业状态的点
     * </p>
     */
    public SplitRoadParams setCheckWorkingStatus(Boolean checkWorkingStatus) {
        this.checkWorkingStatus = checkWorkingStatus;
        return this;
    }

    public Double getPositiveBuffer() {
        return positiveBuffer;
    }

    /**
     * 设置正缓冲，单位（米）
     *
     * @param positiveBuffer
     * @return
     */
    public SplitRoadParams setPositiveBuffer(Double positiveBuffer) {
        this.positiveBuffer = positiveBuffer;
        return this;
    }

    public Double getNegativeBuffer() {
        return negativeBuffer;
    }

    /**
     * 设置负缓冲，单位（米）
     *
     * @param negativeBuffer
     * @return
     */
    public SplitRoadParams setNegativeBuffer(Double negativeBuffer) {
        this.negativeBuffer = negativeBuffer;
        return this;
    }

    public int getAlgorithmIndex() {
        return algorithmIndex;
    }

    public void setAlgorithmIndex(int algorithmIndex) {
        this.algorithmIndex = algorithmIndex;
    }
}