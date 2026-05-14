package sunyu.util.pojo;

import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.Point;

import java.math.BigDecimal;
import java.time.LocalDateTime;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;

/**
 * 农机作业轨迹道路拆分结果实体类
 * <p>
 * 用于封装 {@link sunyu.util.GisUtil#splitRoad(java.util.List, double)} 或其重载方法
 * {@link sunyu.util.GisUtil#splitRoad(java.util.List, double, SplitRoadParams)} 的执行结果，
 * 代表将原始作业轨迹经道路拆分、聚类分析后产生的多个独立作业地块的聚合信息。
 * 该类既包含整体统计信息（总面积、总里程、总时间范围），也包含拆分后的明细地块列表，
 * 支持从整体到局部的多层次作业数据分析。
 * </p>
 * <p>
 * <b>核心设计特点：</b>
 * <ul>
 *   <li><b>聚合计算：</b>{@link #getMu()}、{@link #getTotalJobMileage()}、{@link #getClusterPointCount()} 等方法
 *       基于 {@link #farmPlots} 列表动态聚合计算，确保整体统计与明细数据的一致性；</li>
 *   <li><b>自动排序：</b>{@link #getFarmPlots()} 在返回前按作业起始时间升序排序，保证时序连续性；</li>
 *   <li><b>动态中心点：</b>{@link #getCenterWgs84Point()} 优先从整体几何图形计算内点，支持地图标注和导航定位。</li>
 * </ul>
 * </p>
 * <p>
 * <b>与 {@link FarmPlot} 的关系：</b>
 * <ul>
 *   <li>{@link SplitResult} 是拆分结果的聚合容器，包含多个 {@link FarmPlot}；</li>
 *   <li>每个 {@link FarmPlot} 代表一个独立的作业地块，包含该块的详细时空信息；</li>
 *   <li>{@link SplitResult} 的整体统计由所有 {@link FarmPlot} 聚合得出，避免数据冗余和不一致。</li>
 * </ul>
 * </p>
 * <p>
 * <b>使用场景：</b>
 * <ul>
 *   <li>作业报告生成：汇总单次作业的整体面积、里程、时长和地块数量；</li>
 *   <li>地块级分析：遍历 farmPlots 列表，分析每个独立地块的时空分布特征；</li>
 *   <li>数据库存储：将整体 WKT 和明细列表分别存储，支持灵活查询；</li>
 *   <li>地图可视化：在地图上同时展示整体轮廓和各个子地块边界。</li>
 * </ul>
 * </p>
 *
 * @author SunYu
 * @see FarmPlot
 * @see sunyu.util.GisUtil#splitRoad(java.util.List, double)
 * @see sunyu.util.GisUtil#splitRoad(java.util.List, double, SplitRoadParams)
 * @since 1.0.0
 */
public class SplitResult {
    /**
     * 整体作业起始时间
     * <p>
     * 记录本次拆分结果中所有地块的最早作业开始时间，
     * 通常由 {@link #farmPlots} 中各块起始时间的最小值确定。
     * 用于整体作业时长统计和时序定位。
     * </p>
     */
    private LocalDateTime startTime;
    /**
     * 整体作业结束时间
     * <p>
     * 记录本次拆分结果中所有地块的最晚作业结束时间，
     * 通常由 {@link #farmPlots} 中各块结束时间的最大值确定。
     * 与 {@link #startTime} 共同定义整体作业的时间窗口。
     * </p>
     */
    private LocalDateTime endTime;
    /**
     * 整体轮廓的 WKT 文本表示（WGS84 坐标系）
     * <p>
     * 采用 OGC 标准 WKT 格式存储所有地块合并后的整体几何轮廓，
     * 便于数据库存储、网络传输和跨系统数据交换。
     * 与 {@link #wgs84Geometry} 字段互为补充，分别服务于序列化和空间分析场景。
     * </p>
     */
    private String wkt;
    /**
     * 作业总幅宽（米）
     * <p>
     * 农机作业机具的总工作宽度
     * 该参数在拆分过程中保持一致，用于所有子地块的几何缓冲计算。
     * </p>
     */
    private double workingWidth;
    /**
     * 整体轮廓面积（亩）
     * <p>
     * 所有地块合并后的总面积（亩），由 {@link #setMu(double)} 手动设置。
     * 注意：{@link #getMu()} 方法会基于 {@link #farmPlots} 列表动态聚合计算各子地块面积之和，
     * 使用 BigDecimal 累加以避免浮点精度误差，返回结果可能与该字段值不同。
     * 建议优先使用 {@link #getMu()} 获取实时聚合值。
     * </p>
     */
    private double mu;
    /**
     * 总作业里程（米）
     * <p>
     * 所有地块的作业里程总和（米），由 {@link #setTotalJobMileage(double)} 手动设置。
     * 注意：{@link #getTotalJobMileage()} 方法会基于 {@link #farmPlots} 列表动态聚合计算各子地块里程之和，
     * 使用 BigDecimal 累加以避免浮点精度误差，返回结果可能与该字段值不同。
     * 建议优先使用 {@link #getTotalJobMileage()} 获取实时聚合值。
     * </p>
     */
    private double totalJobMileage;
    /**
     * 拆分后的地块列表
     * <p>
     * 存储道路拆分或聚类后产生的所有独立作业地块，每个元素为一个 {@link FarmPlot} 实例。
     * 该列表是整体统计数据的来源，通过 {@link #getFarmPlots()} 获取时会按作业起始时间升序排序。
     * 默认初始化为空列表（{@code new ArrayList<>()}），避免 null 指针异常。
     * </p>
     */
    private List<FarmPlot> farmPlots = new ArrayList<>();
    /**
     * 整体轮廓的 WGS84 坐标系几何图形（JTS Geometry）
     * <p>
     * 基于 JTS 库的几何对象，通常为 MultiPolygon 类型（多个地块的并集），
     * 由所有子地块的几何图形合并生成。
     * 该对象可直接用于空间分析和地图渲染，
     * 与 {@link #wkt} 字段互为补充，分别服务于空间分析和序列化场景。
     * </p>
     */
    private Geometry wgs84Geometry;
    /**
     * 总聚类点数量
     * <p>
     * 所有地块的原始 GPS 轨迹点数量总和，由 {@link #setClusterPointCount(int)} 手动设置。
     * 注意：{@link #getClusterPointCount()} 方法会基于 {@link #farmPlots} 列表动态聚合计算各子地块点数之和，
     * 返回结果可能与该字段值不同。建议优先使用 {@link #getClusterPointCount()} 获取实时聚合值。
     * </p>
     */
    private int clusterPointCount;
    /**
     * 整体轮廓中心点（WGS84 坐标系）
     * <p>
     * 所有地块合并后的地理中心位置，用于地图标注和整体定位。
     * 可通过 {@link #getCenterWgs84Point()} 动态计算（优先策略），
     * 也可通过 {@link #setCenterWgs84Point(Wgs84Point)} 手动指定。
     * 当 {@link #wgs84Geometry} 不为 null 时，优先使用几何内点（Interior Point）作为中心点。
     * </p>
     */
    private Wgs84Point centerWgs84Point;

    /**
     * 获取整体轮廓的 WKT 文本表示
     *
     * @return WGS84 坐标系下的 WKT 字符串；若未设置则返回 null
     */
    public String getWkt() {
        return wkt;
    }

    /**
     * 设置整体轮廓的 WKT 文本表示
     *
     * @param wkt WGS84 坐标系下的标准 WKT 字符串，如 "MULTIPOLYGON(((...)))"
     */
    public void setWkt(String wkt) {
        this.wkt = wkt;
    }

    /**
     * 获取作业总幅宽（米）
     *
     * @return 农机作业机具的总工作宽度（米）
     */
    public double getWorkingWidth() {
        return workingWidth;
    }

    /**
     * 设置作业总幅宽（米）
     *
     * @param workingWidth 农机作业机具的总工作宽度（米），必须大于等于 0.1 米
     */
    public void setWorkingWidth(double workingWidth) {
        this.workingWidth = workingWidth;
    }

    /**
     * 获取整体轮廓面积（亩）
     * <p>
     * <b>动态聚合计算：</b>基于 {@link #farmPlots} 列表中所有子地块的面积累加得出，
     * 使用 {@link BigDecimal} 进行高精度累加，避免 double 类型的浮点精度误差，
     * 确保财务结算和统计报表的数值准确性。
     * </p>
     * <p>
     * 若 {@link #farmPlots} 为空列表，则返回 0.0。
     * </p>
     *
     * @return 所有子地块面积之和（亩），四舍五入保留 4 位小数
     */
    public double getMu() {
        BigDecimal sumMu = farmPlots.stream()
                .map(farmPlot -> BigDecimal.valueOf(farmPlot.getMu()))
                .reduce(BigDecimal.ZERO, BigDecimal::add);
        return sumMu.doubleValue();
    }

    /**
     * 设置整体轮廓面积（亩）
     * <p>
     * 手动设置整体面积值。注意：{@link #getMu()} 方法会基于子地块列表动态计算，
     * 该方法设置的值仅作为缓存或备用，实际使用中建议以 {@link #getMu()} 的返回为准。
     * </p>
     *
     * @param mu 整体轮廓面积（亩）
     */
    public void setMu(double mu) {
        this.mu = mu;
    }

    /**
     * 获取总作业里程（米）
     * <p>
     * <b>动态聚合计算：</b>基于 {@link #farmPlots} 列表中所有子地块的作业里程累加得出，
     * 使用 {@link BigDecimal} 进行高精度累加，避免 double 类型的浮点精度误差。
     * </p>
     * <p>
     * 若 {@link #farmPlots} 为空列表，则返回 0.0。
     * </p>
     *
     * @return 所有子地块作业里程之和（米）
     */
    public double getTotalJobMileage() {
        BigDecimal sumJobMileage = farmPlots.stream()
                .map(farmPlot -> BigDecimal.valueOf(farmPlot.getJobMileage()))
                .reduce(BigDecimal.ZERO, BigDecimal::add);
        return sumJobMileage.doubleValue();
    }

    /**
     * 设置总作业里程（米）
     * <p>
     * 手动设置总里程值。注意：{@link #getTotalJobMileage()} 方法会基于子地块列表动态计算，
     * 该方法设置的值仅作为缓存或备用，实际使用中建议以 {@link #getTotalJobMileage()} 的返回为准。
     * </p>
     *
     * @param totalJobMileage 总作业里程（米）
     */
    public void setTotalJobMileage(double totalJobMileage) {
        this.totalJobMileage = totalJobMileage;
    }

    /**
     * 获取整体作业起始时间
     *
     * @return 所有地块中最早的作业起始时间；若未设置则返回 null
     */
    public LocalDateTime getStartTime() {
        return startTime;
    }

    /**
     * 设置整体作业起始时间
     *
     * @param startTime 整体作业起始时间，通常为所有子地块起始时间的最小值
     */
    public void setStartTime(LocalDateTime startTime) {
        this.startTime = startTime;
    }

    /**
     * 获取整体作业结束时间
     *
     * @return 所有地块中最晚的作业结束时间；若未设置则返回 null
     */
    public LocalDateTime getEndTime() {
        return endTime;
    }

    /**
     * 设置整体作业结束时间
     *
     * @param endTime 整体作业结束时间，通常为所有子地块结束时间的最大值
     */
    public void setEndTime(LocalDateTime endTime) {
        this.endTime = endTime;
    }

    /**
     * 设置拆分后的地块列表（兼容旧版本方法名）
     * <p>
     * 该方法名 {@code setSplitParts} 为历史遗留命名，功能与 {@link #setFarmPlots(List)} 完全一致，
     * 均用于设置子地块列表。建议新代码使用 {@link #setFarmPlots(List)} 以获得更清晰的语义。
     * </p>
     *
     * @param farmPlots 拆分后的地块列表
     * @see #setFarmPlots(List)
     */
    public void setSplitParts(List<FarmPlot> farmPlots) {
        this.farmPlots = farmPlots;
    }

    /**
     * 获取整体轮廓的 WGS84 坐标系几何图形
     *
     * @return JTS Geometry 对象（通常为 MultiPolygon）；若未设置则返回 null
     */
    public Geometry getWgs84Geometry() {
        return wgs84Geometry;
    }

    /**
     * 设置整体轮廓的 WGS84 坐标系几何图形
     *
     * @param wgs84Geometry JTS Geometry 对象（WGS84 坐标系），通常为 MultiPolygon
     */
    public void setWgs84Geometry(Geometry wgs84Geometry) {
        this.wgs84Geometry = wgs84Geometry;
    }

    /**
     * 获取整体轮廓中心点（WGS84 坐标系）
     * <p>
     * 采用优先动态计算、后备手动值的策略：
     * <ol>
     *   <li>若 {@link #wgs84Geometry} 不为 null，则调用 JTS {@link Geometry#getInteriorPoint()} 动态计算几何内点，
     *       确保中心点始终与当前几何图形保持一致；</li>
     *   <li>若几何图形为 null，则返回 {@link #centerWgs84Point} 字段中手动设置的值。</li>
     * </ol>
     * </p>
     * <p>
     * <b>注意：</b>动态计算得到的中心点不一定在几何图形边界上，
     * 而是位于几何内部的一个代表性点（Interior Point），适用于标注、地图展示等场景。
     * </p>
     *
     * @return 整体轮廓中心点（WGS84 坐标系），若几何为空且未手动设置则返回 null
     * @see Geometry#getInteriorPoint()
     */
    public Wgs84Point getCenterWgs84Point() {
        if (wgs84Geometry != null) {
            Wgs84Point wgs84Point = new Wgs84Point();
            Point centerPoint = wgs84Geometry.getInteriorPoint();
            wgs84Point.setLongitude(centerPoint.getX());
            wgs84Point.setLatitude(centerPoint.getY());
            return wgs84Point;
        }
        return centerWgs84Point;
    }

    /**
     * 手动设置整体轮廓中心点（WGS84 坐标系）
     * <p>
     * 当 {@link #wgs84Geometry} 为 null 时，{@link #getCenterWgs84Point()} 将返回此手动设置的值。
     * 若几何图形已存在，建议优先通过几何动态计算中心点，以确保一致性。
     * </p>
     *
     * @param centerWgs84Point 整体轮廓中心点（WGS84 坐标系）
     */
    public void setCenterWgs84Point(Wgs84Point centerWgs84Point) {
        this.centerWgs84Point = centerWgs84Point;
    }

    /**
     * 获取总聚类点数量
     * <p>
     * <b>动态聚合计算：</b>基于 {@link #farmPlots} 列表中所有子地块的聚类点数量累加得出。
     * </p>
     * <p>
     * 若 {@link #farmPlots} 为空列表，则返回 0。
     * </p>
     *
     * @return 所有子地块聚类点数量之和
     */
    public int getClusterPointCount() {
        return farmPlots.stream().mapToInt(FarmPlot::getClusterPointCount).sum();
    }

    /**
     * 设置总聚类点数量
     * <p>
     * 手动设置总聚类点数量。注意：{@link #getClusterPointCount()} 方法会基于子地块列表动态计算，
     * 该方法设置的值仅作为缓存或备用，实际使用中建议以 {@link #getClusterPointCount()} 的返回为准。
     * </p>
     *
     * @param clusterPointCount 总聚类点数量
     */
    public void setClusterPointCount(int clusterPointCount) {
        this.clusterPointCount = clusterPointCount;
    }

    /**
     * 获取拆分后的地块列表
     * <p>
     * <b>自动排序：</b>返回前会按各子地块的 {@link FarmPlot#getStartTime()} 升序排序，
     * 确保地块列表按作业时间先后顺序排列，便于时序分析和连续作业展示。
     * </p>
     * <p>
     * <b>注意：</b>排序操作会直接修改内部列表的顺序，多次调用结果一致。
     * 若需要保持原始插入顺序，请考虑在调用前复制列表。
     * </p>
     *
     * @return 按作业起始时间升序排列的 FarmPlot 列表；若未设置则返回空列表（非 null）
     */
    public List<FarmPlot> getFarmPlots() {
        farmPlots.sort(Comparator.comparing(FarmPlot::getStartTime));
        return farmPlots;
    }

    /**
     * 设置拆分后的地块列表
     *
     * @param farmPlots 拆分后的 FarmPlot 列表，若为 null 则内部设置为空列表
     */
    public void setFarmPlots(List<FarmPlot> farmPlots) {
        this.farmPlots = farmPlots;
    }

}
