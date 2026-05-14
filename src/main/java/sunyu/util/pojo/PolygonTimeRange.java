package sunyu.util.pojo;

import java.time.LocalDateTime;

/**
 * 多边形（地块）时间范围实体类
 * <p>
 * 用于封装单个多边形（地块）的时间属性，记录该多边形在作业时间轴上的起始和结束时刻。
 * 该类是一个轻量级的不可变（设计意图）数据结构，主要用于时间交叉检测、时序排序和
 * 多边形作业时段的统计分析。
 * </p>
 * <p>
 * <b>核心用途：</b>
 * <ul>
 *   <li>时间重叠检测：判断多个多边形之间是否存在作业时间交叉，为道路拆分和地块合并提供时序依据；</li>
 *   <li>时序排序：按起始时间对多边形进行排序，确保处理顺序的时间连续性；</li>
 *   <li>作业时长统计：计算单个多边形的作业持续时间（end - start），评估作业效率；</li>
 *   <li>时间窗口拆分：当检测到时间重叠时，作为拆分算法的输入参数，指导时间切分策略。</li>
 * </ul>
 * </p>
 * <p>
 * <b>与 {@link TimeRange} 的区别：</b>
 * <ul>
 *   <li>{@link PolygonTimeRange}：关联到具体的多边形索引（polygonIndex），用于多边形级别的时间管理；</li>
 *   <li>{@link TimeRange}：不关联多边形索引，可附带轨迹点集合，用于通用的时间分段和轨迹聚类。</li>
 * </ul>
 * </p>
 * <p>
 * <b>设计说明：</b>
 * 该类仅提供 getter 方法，未提供 setter 方法，设计上倾向于不可变对象，
 * 确保时间范围一旦创建就不会被意外修改，增强多线程环境下的安全性。
 * 若需要修改时间范围，应创建新的实例替换旧实例。
 * </p>
 *
 * @author SunYu
 * @see TimeRange
 * @see sunyu.util.GisUtil#hasTimeOverlap(java.util.Map)
 * @since 1.0.0
 */
public class PolygonTimeRange {
    /**
     * 多边形索引（地块序号）
     * <p>
     * 用于标识该时间范围所属的多边形，与聚类结果或拆分结果中的多边形序号一一对应。
     * 通过该索引可以将时间信息与空间几何信息关联起来，实现时空一体化分析。
     * </p>
     * <p>
     * <b>取值说明：</b>≥ 0 的整数，从 0 开始递增，与 {@link FarmPlotGeometryInfo} 中的键值对应。
     * </p>
     */
    private Integer polygonIndex;
    /**
     * 多边形作业起始时间
     * <p>
     * 记录该多边形（地块）作业开始的最早时间点，
     * 由构成该多边形的所有轨迹点中 GPS 时间的最小值确定。
     * 与 {@link #end} 共同定义该多边形的完整作业时间窗口。
     * </p>
     */
    private LocalDateTime start;
    /**
     * 多边形作业结束时间
     * <p>
     * 记录该多边形（地块）作业结束的最晚时间点，
     * 由构成该多边形的所有轨迹点中 GPS 时间的最大值确定。
     * 与 {@link #start} 共同定义该多边形的完整作业时间窗口。
     * </p>
     * <p>
     * <b>时序约束：</b>在正常情况下，{@code end} 应晚于或等于 {@code start}，
     * 即 {@code end.isAfter(start) || end.isEqual(start)}。
     * 若出现 {@code end} 早于 {@code start} 的情况，表明数据异常，需要校验。
     * </p>
     */
    private LocalDateTime end;

    /**
     * 全参构造方法
     * <p>
     * 创建一个包含多边形索引和完整时间范围的实例。
     * 构造后对象的字段值不可通过 setter 修改，确保时间范围的不可变性。
     * </p>
     *
     * @param polygonIndex 多边形索引（地块序号），必须大于等于 0
     * @param start        多边形作业起始时间，不可为 null
     * @param end          多边形作业结束时间，不可为 null，且应晚于或等于 start
     */
    public PolygonTimeRange(Integer polygonIndex, LocalDateTime start, LocalDateTime end) {
        this.polygonIndex = polygonIndex;
        this.start = start;
        this.end = end;
    }

    /**
     * 获取多边形索引（地块序号）
     *
     * @return 多边形索引（≥ 0 的整数）
     */
    public Integer getPolygonIndex() {
        return polygonIndex;
    }

    /**
     * 获取多边形作业起始时间
     *
     * @return 作业起始时间（LocalDateTime）
     */
    public LocalDateTime getStart() {
        return start;
    }

    /**
     * 获取多边形作业结束时间
     *
     * @return 作业结束时间（LocalDateTime）
     */
    public LocalDateTime getEnd() {
        return end;
    }
}
