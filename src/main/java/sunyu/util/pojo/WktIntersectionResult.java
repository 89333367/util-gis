package sunyu.util.pojo;

/**
 * WKT 几何相交结果实体类
 * <p>
 * 用于封装两个或多个几何图形（通常为地块多边形）空间相交（Intersection）运算后的结果，
 * 包括相交区域的几何形状（WKT 文本）和相交面积（亩）。
 * 该类是空间叠加分析（Overlay Analysis）的基础结果载体，广泛应用于地块重叠检测、
 * 作业区域交叉验证、相邻地块边界冲突分析等场景。
 * </p>
 * <p>
 * <b>核心用途：</b>
 * <ul>
 *   <li>地块重叠检测：判断两个作业地块是否存在空间重叠，并量化重叠面积；</li>
 *   <li>作业区域交叉验证：验证同一农机在不同时间作业的地块是否有交叉，识别重复作业或越界作业；</li>
 *   <li>相邻地块边界分析：计算相邻地块的共享边界长度和重叠面积，辅助边界纠纷处理；</li>
 *   <li>空间查询结果封装：将 JTS 几何相交运算结果转换为可序列化的 WKT + 面积形式，便于存储和传输。</li>
 * </ul>
 * </p>
 * <p>
 * <b>坐标系说明：</b>
 * WKT 文本采用 WGS84 坐标系（EPSG:4326），与系统中其他几何数据保持一致。
 * 面积单位为亩，由 WGS84 几何图形经球面面积计算得出，适用于农业场景的面积统计。
 * </p>
 * <p>
 * <b>与 JTS 的关系：</b>
 * 该类通常由 org.locationtech.jts.geom.Geometry#intersection(Geometry) 运算结果转换而来，
 * 将 JTS Geometry 对象转换为 WKT 文本便于持久化，同时计算并记录相交面积供业务分析使用。
 * </p>
 *
 * @author SunYu
 * @see org.locationtech.jts.geom.Geometry#intersection(org.locationtech.jts.geom.Geometry)
 * @see sunyu.util.GisUtil
 * @since 1.0.0
 */
public class WktIntersectionResult {
    /**
     * 相交区域的几何形状 WKT 文本（WGS84 坐标系）
     * <p>
     * 采用 OGC 标准 WKT 格式存储两个几何图形相交后的结果形状，
     * 几何类型可能为 Polygon（多边形相交）、MultiPolygon（多个多边形相交）、
     * LineString（边界相交）或 Point（顶点相交），具体取决于输入几何的空间关系。
     * </p>
     * <p>
     * <b>典型 WKT 示例：</b>
     * <pre>{@code
     * POLYGON((116.3974 39.9093, 116.3975 39.9093, 116.3975 39.9094, 116.3974 39.9094, 116.3974 39.9093))
     * }</pre>
     * </p>
     * <p>
     * <b>空相交：</b>若两个几何图形不相交（Disjoint），则该字段可能为 null 或空字符串，
     * 此时 {@link #mu} 应为 0.0。
     * </p>
     */
    private String wkt;
    /**
     * 相交区域的面积（亩）
     * <p>
     * 记录两个几何图形相交部分的面积大小，单位为亩。
     * 该面积由 WGS84 坐标系下的几何图形经球面面积计算得出，
     * 反映了实际地理空间中的重叠范围大小。
     * </p>
     * <p>
     * <b>取值说明：</b>
     * <ul>
     *   <li>{@code > 0}：两个几何图形存在有效面积重叠，重叠面积为该值（亩）；</li>
     *   <li>{@code = 0}：两个几何图形不相交或仅边界/顶点接触，无有效面积重叠；</li>
     *   <li>{@code < 0}：异常情况，通常由几何图形自相交或拓扑错误导致，需要校验输入数据。</li>
     * </ul>
     * </p>
     * <p>
     * <b>业务意义：</b>相交面积是判断地块重叠严重程度的核心指标，
     * 面积越大表示重叠越严重，可能需要人工介入处理或调整作业边界。
     * </p>
     */
    private double mu;

    /**
     * 默认无参构造方法
     * <p>
     * 创建一个空的 WktIntersectionResult 实例，wkt 为 null，mu 为 0.0。
     * 通常用于反射创建、序列化反序列化或后续通过 setter 方法填充数据。
     * </p>
     */
    public WktIntersectionResult() {
    }

    /**
     * 全参构造方法
     * <p>
     * 创建一个包含完整相交结果信息的 WktIntersectionResult 实例。
     * 适用于空间相交运算完成后，直接将结果封装为对象的场景。
     * </p>
     *
     * @param wkt 相交区域的 WKT 文本（WGS84 坐标系），若不相交可为 null 或空字符串
     * @param mu  相交区域的面积（亩），不相交时应为 0.0
     */
    public WktIntersectionResult(String wkt, double mu) {
        this.wkt = wkt;
        this.mu = mu;
    }

    /**
     * 获取相交区域的几何形状 WKT 文本
     *
     * @return 相交区域的 WKT 字符串（WGS84 坐标系）；若未设置或不相交则返回 null
     */
    public String getWkt() {
        return wkt;
    }

    /**
     * 设置相交区域的几何形状 WKT 文本
     *
     * @param wkt 相交区域的 WKT 字符串（WGS84 坐标系），如 "POLYGON(((...)))"；
     *            若不相交可设置为 null 或空字符串
     */
    public void setWkt(String wkt) {
        this.wkt = wkt;
    }

    /**
     * 获取相交区域的面积（亩）
     *
     * @return 相交面积（亩）；若未设置则返回 0.0
     */
    public double getMu() {
        return mu;
    }

    /**
     * 设置相交区域的面积（亩）
     *
     * @param mu 相交面积（亩），必须大于等于 0；不相交时应设置为 0.0
     */
    public void setMu(double mu) {
        this.mu = mu;
    }

}
