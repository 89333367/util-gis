package sunyu.util.pojo;

import java.util.List;

/**
 * 时间窗口实体类
 * <p>
 * 用于封装按时间间隔特征分割后的轨迹段，将具有相似时间间隔特征的连续轨迹点归为一个窗口。
 * 该类是时间序列分析和轨迹分段的核心数据结构，通过识别轨迹点之间的时间间隔突变，
 * 将长时段的连续轨迹切分为多个具有内部一致性的子段，便于后续的独立聚类和作业分析。
 * </p>
 * <p>
 * <b>核心设计思想：</b>
 * 农机作业过程中，相邻轨迹点的时间间隔通常较为稳定（如每秒一个点或每 5 秒一个点），
 * 但当农机在田间转移、道路行驶或长时间停车时，时间间隔会出现显著增大。
 * 通过检测这种时间间隔的突变，可以将轨迹自动切分为多个时间窗口，
 * 每个窗口代表一段连续的作业或行驶过程。
 * </p>
 * <p>
 * <b>与 {@link TimeRange} 的区别：</b>
 * <ul>
 *   <li>{@link TimeWindow}：以时间间隔特征（{@link #interval}）为核心划分依据，
 *       携带 WGS84 原始点位，用于轨迹的初步时间分段；</li>
 *   <li>{@link TimeRange}：以绝对时间起止（start/end）为核心，
 *       携带高斯投影点位，用于聚类后的时间范围管理和重叠检测。</li>
 * </ul>
 * </p>
 * <p>
 * <b>典型使用场景：</b>
 * <ul>
 *   <li>轨迹预处理：将原始 GPS 轨迹按时间间隔特征切分为多个窗口，过滤掉停车、怠速等异常段；</li>
 *   <li>作业段识别：识别出时间间隔稳定的连续段作为候选作业段，排除转移和停留段；</li>
 *   <li>分段聚类：对每个时间窗口内的轨迹点独立进行 DBSCAN 空间聚类，避免跨段干扰；</li>
 *   <li>时序分析：统计不同时间窗口的间隔特征，识别作业模式（匀速作业、变速转移、静止停车）。</li>
 * </ul>
 * </p>
 *
 * @author SunYu
 * @see TimeRange
 * @see Wgs84Point
 * @since 1.0.0
 */
public class TimeWindow {
    /**
     * 时间窗口内相邻轨迹点的典型时间间隔（秒）
     * <p>
     * 代表该窗口内所有相邻点位之间时间间隔的特征值（如中位数、均值或众数），
     * 用于描述该轨迹段的时间采样密度和连续性特征。
     * </p>
     * <p>
     * <b>取值含义：</b>
     * <ul>
     *   <li>小值（如 1~10 秒）：表示高频采样，通常为正常作业或行驶状态；</li>
     *   <li>大值（如 60 秒以上）：表示低频采样或时间间隔突变，可能为田间转移、道路行驶或停车；</li>
     *   <li>极大值（如 300 秒以上）：通常表示长时间停车或设备关机，可作为分段边界。</li>
     * </ul>
     * </p>
     * <p>
     * <b>生成方式：</b>由时间窗口分割算法根据窗口内相邻点位的 GPS 时间差统计计算得出，
     * 作为该窗口的时间特征标签。
     * </p>
     */
    private long interval;
    /**
     * 该时间窗口内的 WGS84 坐标轨迹点集合
     * <p>
     * 存储属于该时间窗口的所有原始 GPS 轨迹点，按 GPS 时间升序排列。
     * 这些点位具有相似的时间间隔特征，代表一段连续的作业或行驶过程。
     * </p>
     * <p>
     * <b>数据特点：</b>
     * <ul>
     *   <li>所有点位的 GPS 时间连续，相邻点时间差与 {@link #interval} 特征一致；</li>
     *   <li>坐标系为 WGS84（EPSG:4326），可直接用于地图展示和地理编码；</li>
     *   <li>点位数量反映该窗口的持续时间，可用于判断是否为有效作业段。</li>
     * </ul>
     * </p>
     */
    private List<Wgs84Point> points;

    /**
     * 默认无参构造方法
     * <p>
     * 创建一个空的 TimeWindow 实例，interval 初始化为 0，points 初始化为 null。
     * 通常用于反射创建、序列化反序列化或后续通过 setter 方法填充数据。
     * </p>
     */
    public TimeWindow() {
    }

    /**
     * 全参构造方法
     * <p>
     * 创建一个包含完整时间间隔特征和轨迹点集合的 TimeWindow 实例。
     * 适用于时间窗口分割算法完成后，将分割结果封装为对象的场景。
     * </p>
     *
     * @param interval 时间窗口内相邻轨迹点的典型时间间隔（秒），必须大于等于 0
     * @param points   该窗口内的 WGS84 坐标轨迹点集合，按 GPS 时间升序排列；可为 null
     */
    public TimeWindow(long interval, List<Wgs84Point> points) {
        this.interval = interval;
        this.points = points;
    }

    /**
     * 获取时间窗口内相邻轨迹点的典型时间间隔（秒）
     *
     * @return 时间间隔（秒）；若未设置则返回 0
     */
    public long getInterval() {
        return interval;
    }

    /**
     * 设置时间窗口内相邻轨迹点的典型时间间隔（秒）
     *
     * @param interval 时间间隔（秒），必须大于等于 0
     */
    public void setInterval(long interval) {
        this.interval = interval;
    }

    /**
     * 获取该时间窗口内的 WGS84 坐标轨迹点集合
     *
     * @return WGS84 轨迹点列表，按 GPS 时间升序排列；若未设置则返回 null
     */
    public List<Wgs84Point> getPoints() {
        return points;
    }

    /**
     * 设置该时间窗口内的 WGS84 坐标轨迹点集合
     *
     * @param points WGS84 轨迹点列表，按 GPS 时间升序排列；可为 null
     */
    public void setPoints(List<Wgs84Point> points) {
        this.points = points;
    }
}
