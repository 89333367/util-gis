package sunyu.util.pojo;

import java.time.LocalDateTime;
import java.util.List;

/**
 * 通用时间范围实体类
 * <p>
 * 用于封装一个连续的时间区间（起始时间 → 结束时间），并可选择性地关联该时间范围内的高斯投影轨迹点集合。
 * 该类是农机作业数据分析中的基础数据结构，广泛应用于时间窗口划分、轨迹分段、聚类时序分析和
 * 作业时段统计等场景。
 * </p>
 * <p>
 * <b>与 {@link PolygonTimeRange} 的区别：</b>
 * <ul>
 *   <li>{@link TimeRange}：不关联具体多边形索引，可携带轨迹点集合（{@link #gaussPoints}），
 *       用于通用的时间分段、轨迹聚类和时序分析；</li>
 *   <li>{@link PolygonTimeRange}：关联到具体的多边形索引（polygonIndex），不可携带轨迹点，
 *       用于多边形级别的时间管理和重叠检测。</li>
 * </ul>
 * </p>
 * <p>
 * <b>核心使用场景：</b>
 * <ul>
 *   <li>时间窗口拆分：将连续作业轨迹按时间间隙切分为多个不重叠的时间段；</li>
 *   <li>轨迹分段聚类：对每个时间窗口内的轨迹点进行独立的空间聚类分析；</li>
 *   <li>作业时段统计：计算每个时间段的作业时长、平均速度、里程等指标；</li>
 *   <li>时序重叠检测：判断多个时间范围之间是否存在交叉或包含关系。</li>
 * </ul>
 * </p>
 * <p>
 * <b>时序约束：</b>在正常情况下，{@code end} 应晚于或等于 {@code start}，
 * 即 {@code end.isAfter(start) || end.isEqual(start)}。
 * 若出现 {@code end} 早于 {@code start} 的情况，表明数据异常，需要校验。
 * </p>
 *
 * @author SunYu
 * @see PolygonTimeRange
 * @see TimeWindow
 * @see GaussPoint
 * @since 1.0.0
 */
public class TimeRange {
    /**
     * 时间范围起始时刻
     * <p>
     * 记录该时间段的开始时间点，通常由轨迹点集合中 GPS 时间的最小值确定，
     * 或由时间拆分算法根据时间间隙自动计算得出。
     * 与 {@link #end} 共同定义一个闭合的时间区间 [{@code start}, {@code end}]。
     * </p>
     */
    private LocalDateTime start;
    /**
     * 时间范围结束时刻
     * <p>
     * 记录该时间段的结束时间点，通常由轨迹点集合中 GPS 时间的最大值确定，
     * 或由时间拆分算法根据时间间隙自动计算得出。
     * 与 {@link #start} 共同定义一个闭合的时间区间 [{@code start}, {@code end}]。
     * </p>
     * <p>
     * <b>时序约束：</b>正常情况下应满足 {@code end.isAfter(start) || end.isEqual(start)}。
     * </p>
     */
    private LocalDateTime end;
    /**
     * 该时间范围内的轨迹点集合（高斯投影坐标系）
     * <p>
     * 存储属于该时间区间的所有高斯投影轨迹点，按 GPS 时间升序排列。
     * 该字段为可选属性，在仅需要时间范围而不关心具体轨迹点的场景下可为 null。
     * </p>
     * <p>
     * <b>典型使用场景：</b>
     * <ul>
     *   <li>时间窗口拆分后，每个时间范围携带对应的轨迹点子集，用于后续独立聚类；</li>
     *   <li>作业时段分析时，提取该时间段内的所有轨迹点进行速度和里程计算。</li>
     * </ul>
     * </p>
     */
    private List<GaussPoint> gaussPoints;

    /**
     * 全参构造方法（含轨迹点集合）
     * <p>
     * 创建一个包含完整时间信息和关联轨迹点集合的 TimeRange 实例。
     * 适用于时间窗口拆分后需要将轨迹点与时间范围一并封装的场景。
     * </p>
     *
     * @param start       时间范围起始时刻，不可为 null
     * @param end         时间范围结束时刻，不可为 null，应晚于或等于 start
     * @param gaussPoints 该时间范围内的轨迹点集合（高斯投影坐标系），可为 null
     */
    public TimeRange(LocalDateTime start, LocalDateTime end, List<GaussPoint> gaussPoints) {
        this.start = start;
        this.end = end;
        this.gaussPoints = gaussPoints;
    }

    /**
     * 部分参构造方法（不含轨迹点集合）
     * <p>
     * 创建一个仅包含时间起止信息的 TimeRange 实例，轨迹点集合初始化为 null。
     * 适用于仅需时间范围进行重叠检测或排序，而无需关联具体轨迹点的场景。
     * </p>
     *
     * @param start 时间范围起始时刻，不可为 null
     * @param end   时间范围结束时刻，不可为 null，应晚于或等于 start
     */
    public TimeRange(LocalDateTime start, LocalDateTime end) {
        this.start = start;
        this.end = end;
    }

    /**
     * 默认无参构造方法
     * <p>
     * 创建一个空的 TimeRange 实例，所有字段初始化为 null。
     * 通常用于反射创建、序列化反序列化或后续通过 setter 方法填充数据。
     * </p>
     */
    public TimeRange() {
    }

    /**
     * 获取时间范围起始时刻
     *
     * @return 起始时间（LocalDateTime）；若未设置则返回 null
     */
    public LocalDateTime getStart() {
        return start;
    }

    /**
     * 设置时间范围起始时刻
     *
     * @param start 起始时间（LocalDateTime），不可为 null
     */
    public void setStart(LocalDateTime start) {
        this.start = start;
    }

    /**
     * 获取时间范围结束时刻
     *
     * @return 结束时间（LocalDateTime）；若未设置则返回 null
     */
    public LocalDateTime getEnd() {
        return end;
    }

    /**
     * 设置时间范围结束时刻
     *
     * @param end 结束时间（LocalDateTime），不可为 null，应晚于或等于 start
     */
    public void setEnd(LocalDateTime end) {
        this.end = end;
    }

    /**
     * 获取该时间范围内的轨迹点集合
     *
     * @return 高斯投影轨迹点列表；若未设置则返回 null
     */
    public List<GaussPoint> getGaussPoints() {
        return gaussPoints;
    }

    /**
     * 设置该时间范围内的轨迹点集合
     *
     * @param gaussPoints 高斯投影轨迹点列表，按 GPS 时间升序排列；可为 null
     */
    public void setGaussPoints(List<GaussPoint> gaussPoints) {
        this.gaussPoints = gaussPoints;
    }
}
