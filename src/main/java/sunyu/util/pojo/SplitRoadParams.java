package sunyu.util.pojo;

/**
 * 农机作业轨迹道路拆分算法参数配置实体类
 * <p>
 * 用于封装 {@link sunyu.util.GisUtil#splitRoad(java.util.List, double, SplitRoadParams)} 方法的可选参数，
 * 提供对 DBSCAN 空间聚类、几何缓冲、面积过滤、作业状态过滤和拆分策略等算法的精细化控制。
 * 该类采用 Builder 风格设计（setter 返回 {@code this}），支持链式调用，便于参数的快速配置和组合。
 * </p>
 * <p>
 * <b>核心参数分组：</b>
 * <ul>
 *   <li><b>DBSCAN 聚类参数：</b>{@link #dbScanEpsilon}（邻域半径）和 {@link #dbScanMinPoints}（最小点数），
 *       控制空间聚类的密度阈值，决定哪些轨迹点被归为同一地块；</li>
 *   <li><b>几何缓冲参数：</b>{@link #positiveBuffer}（正向缓冲）和 {@link #negativeBuffer}（负向缓冲），
 *       用于调整地块边界的扩张或收缩，消除锯齿和填补缝隙；</li>
 *   <li><b>面积过滤参数：</b>{@link #minReturnMu}（最小返回面积），过滤掉面积过小的无效地块；</li>
 *   <li><b>作业状态参数：</b>{@link #checkWorkingStatus}，控制是否过滤未作业状态的轨迹点；</li>
 *   <li><b>拆分策略参数：</b>{@link #algorithmIndex}，选择空间合并或时间拆分两种处理模式。</li>
 * </ul>
 * </p>
 * <p>
 * <b>使用示例：</b>
 * <pre>{@code
 * SplitRoadParams params = new SplitRoadParams()
 *     .setDbScanEpsilon(11.0)
 *     .setDbScanMinPoints(30)
 *     .setMinReturnMu(0.55)
 *     .setCheckWorkingStatus(true)
 *     .setAlgorithmIndex(1);
 * SplitResult result = gisUtil.splitRoad(points, 3.5, params);
 * }</pre>
 * </p>
 * <p>
 * <b>默认值说明：</b>
 * <ul>
 *   <li>{@link #checkWorkingStatus}：默认 {@code false}，不过滤未作业状态点；</li>
 *   <li>{@link #algorithmIndex}：默认 {@code 1}（时间算法），优先保证时序连续性。</li>
 * </ul>
 * </p>
 *
 * @author SunYu
 * @see sunyu.util.GisUtil#splitRoad(java.util.List, double, SplitRoadParams)
 * @see SplitResult
 * @since 1.0.0
 */
public class SplitRoadParams {
    /**
     * DBSCAN 聚类邻域半径（米）
     * <p>
     * 用于 DBSCAN（Density-Based Spatial Clustering of Applications with Noise）算法的核心参数 ε（epsilon），
     * 定义了一个点周围邻域的半径范围。当两个轨迹点之间的欧几里得距离小于此值时，
     * 它们被视为邻居点（密度可达），属于同一个聚类簇（地块）。
     * </p>
     * <p>
     * <b>取值建议：</b>
     * <ul>
     *   <li>值过小：导致聚类过于细碎，同一地块被拆分为多个小块；</li>
     *   <li>值过大：导致不同地块被错误合并，丢失道路边界信息；</li>
     *   <li>推荐范围：通常根据农机作业幅宽和 GPS 精度设置，一般为 10~30 米。</li>
     * </ul>
     * </p>
     * <p>
     * <b>与 {@link #dbScanMinPoints} 的关系：</b>两者共同决定聚类密度，
     * 只有当邻域内点数 ≥ dbScanMinPoints 时，该点才会被确认为核心点并扩展聚类。
     * </p>
     */
    private Double dbScanEpsilon;
    /**
     * DBSCAN 聚类最小点数
     * <p>
     * 用于 DBSCAN 算法的核心参数 MinPts，定义了一个邻域内被认定为高密度区域所需的最少点数。
     * 当一个点 ε 邻域内的轨迹点数量大于等于此值时，该点被标记为核心点（Core Point），
     * 并以此为中心向外扩展聚类；反之，若邻域内点数不足，则该点被视为噪声点或边界点，
     * 不会被分配到任何聚类中（最终被过滤）。
     * </p>
     * <p>
     * <b>取值建议：</b>
     * <ul>
     *   <li>值过小：噪声点容易被误判为有效聚类，产生大量碎小地块；</li>
     *   <li>值过大：有效稀疏地块可能被误判为噪声，导致真实作业区域丢失；</li>
     *   <li>推荐范围：通常设置为 3~10，结合轨迹点采样频率和作业速度调整。</li>
     * </ul>
     * </p>
     */
    private Integer dbScanMinPoints;
    /**
     * 正向几何缓冲距离（米）
     * <p>
     * 对聚类生成的地块多边形进行向外扩张（缓冲）的距离，用于填补地块边缘的锯齿缝隙、
     * 平滑边界或补偿 GPS 定位误差造成的轮廓缺失。
     * 正值表示向外扩张，负值表示向内收缩（通常配合 {@link #negativeBuffer} 使用）。
     * </p>
     * <p>
     * <b>典型使用场景：</b>
     * <ul>
     *   <li>先正向缓冲（如 +2 米）填补缝隙，再负向缓冲（如 -2 米）恢复原始尺寸，
     *       实现边界平滑而不改变整体面积；</li>
     *   <li>仅正向缓冲以扩大作业区域，覆盖因 GPS 漂移导致的边缘遗漏。</li>
     * </ul>
     * </p>
     */
    private Double positiveBuffer;
    /**
     * 负向几何缓冲距离（米）
     * <p>
     * 对聚类生成的地块多边形进行向内收缩（缓冲）的距离，用于消除地块边缘的毛刺、
     * 剔除因 GPS 噪声产生的冗余突出部分，或抵消正向缓冲带来的面积膨胀。
     * 负值表示向内收缩，通常与 {@link #positiveBuffer} 配合使用以达到边界优化效果。
     * </p>
     * <p>
     * <b>典型使用场景：</b>
     * <ul>
     *   <li>与正向缓冲组合使用（如先 +2 米再 -2 米），实现边界平滑和孔洞填补；</li>
     *   <li>单独使用以剔除地块边缘的细碎突出，使轮廓更加规整。</li>
     * </ul>
     * </p>
     */
    private Double negativeBuffer;
    /**
     * 最小返回地块面积阈值（亩）
     * <p>
     * 用于过滤面积过小的无效地块。在聚类和缓冲操作完成后，
     * 所有生成的地块多边形会按面积进行筛选，面积小于此阈值的地块将被丢弃，
     * 不会出现在最终的 {@link SplitResult#getFarmPlots()} 结果中。
     * </p>
     * <p>
     * <b>取值建议：</b>
     * <ul>
     *   <li>值过小：保留大量碎小地块（如转弯掉头区域），影响数据质量；</li>
     *   <li>值过大：可能误过滤掉真实的小型作业地块；</li>
     *   <li>推荐范围：根据农机类型和作业精度要求设置，一般为 0.1~1.0 亩。</li>
     * </ul>
     * </p>
     */
    private Double minReturnMu;
    /**
     * 是否根据作业状态过滤轨迹点
     * <p>
     * 控制聚类计算前是否先过滤掉未作业状态（非作业、怠速、停车等）的轨迹点。
     * 当设置为 {@code true} 时，只有处于作业状态的点才会参与 DBSCAN 聚类和后续几何计算，
     * 有效排除道路行驶、田间转移等非作业轨迹对地块轮廓的干扰。
     * </p>
     * <p>
     * <b>取值说明：</b>
     * <ul>
     *   <li>{@code true}：启用作业状态过滤，仅保留作业状态点参与聚类，结果更精确但可能丢失部分边界信息；</li>
     *   <li>{@code false}（默认）：不启用过滤，所有轨迹点均参与聚类，结果更完整但可能包含非作业区域。</li>
     * </ul>
     * </p>
     */
    private Boolean checkWorkingStatus = false;

    /**
     * 道路拆分算法策略索引
     * <p>
     * 控制当检测到多个地块之间存在时间重叠时的处理策略，提供两种互斥的算法模式：
     * </p>
     * <p>
     * <b>取值说明：</b>
     * <ul>
     *   <li>{@code 0} — <b>空间算法（空间合并模式）：</b>
     *     先使用 DBSCAN 进行空间聚类，若聚类后的不同地块之间存在时间交叉，
     *     则将时间交叉的地块进行空间合并，生成一个更大的联合地块。
     *     适用于同一田块内多次往返作业、时间重叠但空间连续的场景；</li>
     *   <li>{@code 1}（默认）— <b>时间算法（时间拆分模式）：</b>
     *     先使用 DBSCAN 进行空间聚类，再对每个聚类簇按时间窗口进行二次拆分，
     *     若地块之间存在时间交叉，则在时间交叉的地块内部再次进行时间窗口拆分，
     *     确保最终每个地块具有独立且不重叠的时间范围。
     *     适用于跨田块作业、道路转移导致时间重叠但空间分离的场景。</li>
     * </ul>
     * </p>
     * <p>
     * <b>选择建议：</b>
     * <ul>
     *   <li>若作业轨迹在同一地块内多次往返（如犁地、播种），时间有重叠但空间连续，选择 {@code 0}；</li>
     *   <li>若作业轨迹跨越多个独立地块（如地块 A → 道路 → 地块 B），时间有重叠但空间分离，选择 {@code 1}。</li>
     * </ul>
     * </p>
     */
    private int algorithmIndex = 1;

    /**
     * 默认无参构造方法
     * <p>
     * 创建一个使用默认参数配置的 SplitRoadParams 实例：
     * <ul>
     *   <li>{@link #checkWorkingStatus} = {@code false}</li>
     *   <li>{@link #algorithmIndex} = {@code 1}（时间算法）</li>
     *   <li>其他数值参数均为 {@code null}，表示使用算法内部默认值</li>
     * </ul>
     * </p>
     */
    public SplitRoadParams() {
    }

    /**
     * 获取 DBSCAN 聚类邻域半径（米）
     *
     * @return 邻域半径 ε（米）；若未设置则返回 null，表示使用算法内部默认值
     */
    public Double getDbScanEpsilon() {
        return dbScanEpsilon;
    }

    /**
     * 设置 DBSCAN 聚类邻域半径（米）
     * <p>
     * 设置 DBSCAN 算法的 ε 参数，控制空间聚类的密度半径。
     * 两个轨迹点之间的距离小于此值时被视为同一聚类的邻居。
     * </p>
     *
     * @param dbScanEpsilon 邻域半径（米），必须大于 0，推荐 10~30 米
     * @return 当前实例（{@code this}），支持链式调用
     */
    public SplitRoadParams setDbScanEpsilon(Double dbScanEpsilon) {
        this.dbScanEpsilon = dbScanEpsilon;
        return this;
    }

    /**
     * 获取 DBSCAN 聚类最小点数
     *
     * @return 最小点数 MinPts；若未设置则返回 null，表示使用算法内部默认值
     */
    public Integer getDbScanMinPoints() {
        return dbScanMinPoints;
    }

    /**
     * 设置 DBSCAN 聚类最小点数
     * <p>
     * 设置 DBSCAN 算法的 MinPts 参数，控制形成高密度区域所需的最少点数。
     * 邻域内点数不足此值的点将被视为噪声并过滤。
     * </p>
     *
     * @param dbScanMinPoints 最小点数，必须大于等于 1，推荐 3~10
     * @return 当前实例（{@code this}），支持链式调用
     */
    public SplitRoadParams setDbScanMinPoints(Integer dbScanMinPoints) {
        this.dbScanMinPoints = dbScanMinPoints;
        return this;
    }

    /**
     * 获取最小返回地块面积阈值（亩）
     *
     * @return 最小面积阈值（亩）；若未设置则返回 null，表示不过滤
     */
    public Double getMinReturnMu() {
        return minReturnMu;
    }

    /**
     * 设置最小返回地块面积阈值（亩）
     * <p>
     * 面积小于此阈值的地块将在最终结果中被过滤丢弃，
     * 用于排除因 GPS 噪声或短暂停留产生的碎小无效地块。
     * </p>
     *
     * @param minReturnMu 最小面积阈值（亩），必须大于等于 0，推荐 0.1~1.0 亩
     * @return 当前实例（{@code this}），支持链式调用
     */
    public SplitRoadParams setMinReturnMu(Double minReturnMu) {
        this.minReturnMu = minReturnMu;
        return this;
    }

    /**
     * 获取是否根据作业状态过滤轨迹点
     *
     * @return {@code true} 表示启用作业状态过滤；{@code false}（默认）表示不过滤
     */
    public Boolean getCheckWorkingStatus() {
        return checkWorkingStatus;
    }

    /**
     * 设置是否根据作业状态过滤轨迹点
     * <p>
     * 启用后（{@code true}），只有作业状态的轨迹点会参与聚类和几何计算，
     * 有效排除道路行驶、田间转移等非作业轨迹对地块轮廓的干扰。
     * </p>
     *
     * @param checkWorkingStatus {@code true} 启用过滤，{@code false} 禁用过滤
     * @return 当前实例（{@code this}），支持链式调用
     */
    public SplitRoadParams setCheckWorkingStatus(Boolean checkWorkingStatus) {
        this.checkWorkingStatus = checkWorkingStatus;
        return this;
    }

    /**
     * 获取正向几何缓冲距离（米）
     *
     * @return 正向缓冲距离（米）；若未设置则返回 null
     */
    public Double getPositiveBuffer() {
        return positiveBuffer;
    }

    /**
     * 设置正向几何缓冲距离（米）
     * <p>
     * 对地块多边形进行向外扩张的距离，用于填补边缘缝隙、平滑边界或补偿 GPS 误差。
     * 正值表示向外扩张，通常与 {@link #negativeBuffer} 配合使用。
     * </p>
     *
     * @param positiveBuffer 正向缓冲距离（米），正值表示向外扩张
     * @return 当前实例（{@code this}），支持链式调用
     */
    public SplitRoadParams setPositiveBuffer(Double positiveBuffer) {
        this.positiveBuffer = positiveBuffer;
        return this;
    }

    /**
     * 获取负向几何缓冲距离（米）
     *
     * @return 负向缓冲距离（米）；若未设置则返回 null
     */
    public Double getNegativeBuffer() {
        return negativeBuffer;
    }

    /**
     * 设置负向几何缓冲距离（米）
     * <p>
     * 对地块多边形进行向内收缩的距离，用于消除边缘毛刺、剔除冗余突出部分，
     * 或抵消正向缓冲带来的面积膨胀。负值表示向内收缩，通常与 {@link #positiveBuffer} 配合使用。
     * </p>
     *
     * @param negativeBuffer 负向缓冲距离（米），负值表示向内收缩
     * @return 当前实例（{@code this}），支持链式调用
     */
    public SplitRoadParams setNegativeBuffer(Double negativeBuffer) {
        this.negativeBuffer = negativeBuffer;
        return this;
    }

    /**
     * 获取道路拆分算法策略索引
     *
     * @return 算法索引，{@code 0} 表示空间算法（合并模式），{@code 1}（默认）表示时间算法（拆分模式）
     */
    public int getAlgorithmIndex() {
        return algorithmIndex;
    }

    /**
     * 设置道路拆分算法策略索引
     * <p>
     * 控制时间重叠地块的处理策略：
     * <ul>
     *   <li>{@code 0}：空间算法，时间重叠的地块进行空间合并；</li>
     *   <li>{@code 1}（默认）：时间算法，时间重叠的地块进行时间窗口拆分。</li>
     * </ul>
     * </p>
     *
     * @param algorithmIndex 算法索引，{@code 0} 或 {@code 1}
     * @return 当前实例（{@code this}），支持链式调用
     */
    public SplitRoadParams setAlgorithmIndex(int algorithmIndex) {
        this.algorithmIndex = algorithmIndex;
        return this;
    }
}
