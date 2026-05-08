package sunyu.util.pojo;

import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.Point;

import java.time.LocalDateTime;
import java.util.List;

/**
 * 农机作业地块信息实体类
 * <p>
 * 用于封装农机在田间作业过程中形成的单个作业地块的完整信息，
 * 包括作业时间范围、作业面积、作业里程、几何图形、中心点坐标及关联的轨迹点集合等核心属性。
 * </p>
 * <p>
 * 该类是 {@link SplitResult} 的组成单元，代表道路拆分或聚类后的一个独立作业区域。
 * 几何数据统一基于 WGS84 坐标系（EPSG:4326），面积单位为亩，距离单位为米。
 * </p>
 * <p>
 * <b>核心属性说明：</b>
 * <ul>
 *   <li>时间属性：记录该地块作业的起止时间，用于作业时长统计和时序分析</li>
 *   <li>面积属性：以亩为单位的地块面积，由 WGS84 几何图形经球面面积计算得出</li>
 *   <li>里程属性：农机在该地块内的实际行驶里程（米），基于高斯投影坐标累加计算</li>
 *   <li>几何属性：WKT 文本和 JTS Geometry 对象双重存储，兼顾序列化和空间分析</li>
 *   <li>中心点：通过 {@link #getCenterWgs84Point()} 动态计算几何内点</li>
 *   <li>轨迹点：保留构成该地块的原始高斯投影轨迹点，支持溯源和精细分析</li>
 * </ul>
 * </p>
 *
 * @author SunYu
 * @see SplitResult
 * @see Wgs84Point
 * @see GaussPoint
 * @since 1.0.0
 */
public class FarmPlot {
    /**
     * 地块作业起始时间
     * <p>
     * 记录农机进入该地块开始作业的最早时间点，
     * 由构成该地块的轨迹点集合中 GPS 时间的最小值确定。
     * 用于作业时长统计、时序分析和与其他地块的先后顺序判断。
     * </p>
     */
    private LocalDateTime startTime;
    /**
     * 地块作业结束时间
     * <p>
     * 记录农机离开该地块结束作业的最晚时间点，
     * 由构成该地块的轨迹点集合中 GPS 时间的最大值确定。
     * 与 {@link #startTime} 共同确定该地块的作业时间窗口。
     * </p>
     */
    private LocalDateTime endTime;
    /**
     * 地块面积（亩）
     * <p>
     * 该地块在 WGS84 坐标系下的球面面积，经国家标准换算为亩（1 亩 = 2000/3 平方米）。
     * 由 {@link sunyu.util.GisUtil#calcMu(Geometry)} 方法基于几何图形计算得出，
     * 四舍五入保留 4 位小数，满足财务结算精度要求。
     * </p>
     */
    private double mu;
    /**
     * 作业里程（米）
     * <p>
     * 农机在该地块内实际行驶轨迹的总长度，基于高斯投影坐标系下相邻轨迹点的欧几里得距离累加计算。
     * 区别于地块几何图形的周长，作业里程反映的是实际行驶路径长度，
     * 用于油耗估算、作业效率评估和设备损耗统计。
     * </p>
     */
    private double jobMileage;
    /**
     * 地块几何图形的 WKT 文本表示（WGS84 坐标系）
     * <p>
     * 采用 OGC 标准 WKT（Well-Known Text）格式存储地块的几何轮廓，
     * 典型格式如 "POLYGON((经度1 纬度1, 经度2 纬度2, ...))"。
     * 该字段便于数据库存储、网络传输和跨系统数据交换，
     * 与 {@link #wgs84Geometry} 字段互为补充，分别服务于序列化和空间分析场景。
     * </p>
     */
    private String wkt;
    /**
     * 地块的 WGS84 坐标系几何图形（JTS Geometry）
     * <p>
     * 基于 JTS（Java Topology Suite）库的几何对象，通常为 Polygon 或 MultiPolygon 类型，
     * 由作业轨迹经缓冲、合并、简化等空间操作生成。
     * 该对象可直接用于空间分析（如相交、包含、面积计算等），
     * 与 {@link #wkt} 字段互为补充，分别服务于空间分析和序列化场景。
     * </p>
     */
    private Geometry wgs84Geometry;
    /**
     * 聚类点数量
     * <p>
     * 构成该地块的原始 GPS 轨迹点数量，即 DBSCAN 聚类或道路拆分后归入该地块的点位总数。
     * 该值反映地块的数据密度和轨迹覆盖程度，
     * 用于数据质量评估（如判断是否存在点位稀疏导致的面积计算偏差）。
     * </p>
     */
    private int clusterPointCount;
    /**
     * 作业总幅宽（米）
     * <p>
     * 农机作业机具的总工作宽度，即单次作业覆盖的地面宽度。
     * 该参数直接影响地块几何图形的缓冲半径和面积计算结果，
     * 是连接轨迹线与实际作业区域的关键参数。
     * 例如：播种机幅宽 4 米，则轨迹线两侧各缓冲 2 米形成作业区域。
     * </p>
     */
    private double workingWidth;
    /**
     * 地块中心点（WGS84 坐标系）
     * <p>
     * 该地块的地理中心位置，用于地图标注、导航定位和空间索引。
     * 可通过 {@link #getCenterWgs84Point()} 方法动态计算（优先策略），
     * 也可通过 {@link #setCenterWgs84Point(Wgs84Point)} 手动指定。
     * 当 {@link #wgs84Geometry} 不为 null 时，优先使用几何内点（Interior Point）作为中心点。
     * </p>
     */
    private Wgs84Point centerWgs84Point;

    /**
     * 构成该地块的原始高斯投影轨迹点集合
     * <p>
     * 保留构成该地块的所有原始 GPS 轨迹点（已转换为高斯投影坐标系），
     * 按 GPS 时间先后顺序排列。
     * 该集合支持轨迹溯源、作业里程精确计算、时间窗口分析和轨迹回放等高级功能，
     * 是连接几何结果与原始数据的关键纽带。
     * </p>
     */
    private List<GaussPoint> geometryPoints;

    /**
     * 获取构成该地块的原始高斯投影轨迹点集合
     *
     * @return 按 GPS 时间排序的 GaussPoint 列表；若未设置则返回 null
     */
    public List<GaussPoint> getGeometryPoints() {
        return geometryPoints;
    }

    /**
     * 设置构成该地块的原始高斯投影轨迹点集合
     *
     * @param geometryPoints 按 GPS 时间排序的 GaussPoint 列表
     */
    public void setGeometryPoints(List<GaussPoint> geometryPoints) {
        this.geometryPoints = geometryPoints;
    }

    /**
     * 获取地块作业起始时间
     *
     * @return 作业起始时间；若未设置则返回 null
     */
    public LocalDateTime getStartTime() {
        return startTime;
    }

    /**
     * 设置地块作业起始时间
     *
     * @param startTime 作业起始时间，通常为该地块最早轨迹点的 GPS 时间
     */
    public void setStartTime(LocalDateTime startTime) {
        this.startTime = startTime;
    }

    /**
     * 获取地块作业结束时间
     *
     * @return 作业结束时间；若未设置则返回 null
     */
    public LocalDateTime getEndTime() {
        return endTime;
    }

    /**
     * 设置地块作业结束时间
     *
     * @param endTime 作业结束时间，通常为该地块最晚轨迹点的 GPS 时间
     */
    public void setEndTime(LocalDateTime endTime) {
        this.endTime = endTime;
    }

    /**
     * 获取地块面积（亩）
     *
     * @return 地块面积（亩），四舍五入保留 4 位小数；若未计算则返回 0.0
     */
    public double getMu() {
        return mu;
    }

    /**
     * 设置地块面积（亩）
     *
     * @param mu 地块面积（亩），通常由 {@link sunyu.util.GisUtil#calcMu(Geometry)} 计算得出
     */
    public void setMu(double mu) {
        this.mu = mu;
    }

    /**
     * 获取作业里程（米）
     *
     * @return 农机在该地块内的实际行驶里程（米）；若未计算则返回 0.0
     */
    public double getJobMileage() {
        return jobMileage;
    }

    /**
     * 设置作业里程（米）
     *
     * @param jobMileage 作业里程（米），基于高斯投影坐标累加计算
     */
    public void setJobMileage(double jobMileage) {
        this.jobMileage = jobMileage;
    }

    /**
     * 获取地块几何图形的 WKT 文本表示
     *
     * @return WGS84 坐标系下的 WKT 字符串；若未设置则返回 null
     */
    public String getWkt() {
        return wkt;
    }

    /**
     * 设置地块几何图形的 WKT 文本表示
     *
     * @param wkt WGS84 坐标系下的标准 WKT 字符串，如 "POLYGON((...))"
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
     * 获取聚类点数量
     *
     * @return 构成该地块的原始 GPS 轨迹点数量
     */
    public int getClusterPointCount() {
        return clusterPointCount;
    }

    /**
     * 设置聚类点数量
     *
     * @param clusterPointCount 构成该地块的原始 GPS 轨迹点数量
     */
    public void setClusterPointCount(int clusterPointCount) {
        this.clusterPointCount = clusterPointCount;
    }

    /**
     * 获取地块的 WGS84 坐标系几何图形
     *
     * @return JTS Geometry 对象（通常为 Polygon 或 MultiPolygon）；若未设置则返回 null
     */
    public Geometry getWgs84Geometry() {
        return wgs84Geometry;
    }

    /**
     * 设置地块的 WGS84 坐标系几何图形
     *
     * @param wgs84Geometry JTS Geometry 对象（WGS84 坐标系），通常为 Polygon 或 MultiPolygon
     */
    public void setWgs84Geometry(Geometry wgs84Geometry) {
        this.wgs84Geometry = wgs84Geometry;
    }

    /**
     * 获取地块的中心点坐标（WGS84坐标系）
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
     * @return 地块中心点（WGS84坐标系），若几何为空且未手动设置则返回 null
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
     * 手动设置地块中心点（WGS84 坐标系）
     * <p>
     * 当 {@link #wgs84Geometry} 为 null 时，{@link #getCenterWgs84Point()} 将返回此手动设置的值。
     * 若几何图形已存在，建议优先通过几何动态计算中心点，以确保一致性。
     * </p>
     *
     * @param centerWgs84Point 地块中心点（WGS84 坐标系）
     */
    public void setCenterWgs84Point(Wgs84Point centerWgs84Point) {
        this.centerWgs84Point = centerWgs84Point;
    }

}
