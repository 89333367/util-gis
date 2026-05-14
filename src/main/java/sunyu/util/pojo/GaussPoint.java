package sunyu.util.pojo;

/**
 * 高斯-克吕格投影坐标点实体类
 * <p>
 * 继承自 {@link Wgs84Point}，在保留原始 WGS84 经纬度坐标的基础上，
 * 扩展了高斯-克吕格投影（Gauss-Krüger）平面直角坐标（X, Y），
 * 实现了一个点位同时具有地理坐标系（WGS84）和投影坐标系（高斯）的双重表示。
 * </p>
 * <p>
 * <b>坐标系说明：</b>
 * <ul>
 *   <li>WGS84 坐标：继承自父类的经度（{@link Wgs84Point#longitude}）和纬度（{@link Wgs84Point#latitude}），
 *       用于全球定位、地图展示和跨系统数据交换；</li>
 *   <li>高斯投影坐标：本类扩展的 {@link #gaussX}（北向坐标，单位：米）和 {@link #gaussY}（东向坐标，单位：米），
 *       用于平面距离计算、几何缓冲、面积量算等空间分析操作。</li>
 * </ul>
 * </p>
 * <p>
 * <b>核心用途：</b>
 * <ul>
 *   <li>轨迹点投影转换：将 GPS 采集的 WGS84 经纬度点转换为高斯投影平面坐标，消除地球曲率影响；</li>
 *   <li>距离精确计算：基于平面直角坐标计算欧几里得距离，精度远高于球面距离公式；</li>
 *   <li>几何操作输入：作为 JTS 几何构建、缓冲区生成、多边形合并等空间算法的坐标源；</li>
 *   <li>聚类与拆分：在 DBSCAN 聚类和道路拆分算法中，使用高斯坐标进行密度聚类和空间分段。</li>
 * </ul>
 * </p>
 * <p>
 * <b>重要说明：</b>
 * 高斯投影坐标（gaussX, gaussY）由 {@link sunyu.util.GisUtil#toGaussPointList(java.util.List)} 等工具方法生成，
 * 基于几何中心经度自动计算 6° 分带投影参数（带号、中央经线、假东距）。
 * 同一批轨迹点应使用统一的投影带参数，确保坐标一致性。
 * </p>
 *
 * @author SunYu
 * @see Wgs84Point
 * @see sunyu.util.GisUtil#toGaussPointList(java.util.List)
 * @since 1.0.0
 */
public class GaussPoint extends Wgs84Point {
    /**
     * 高斯投影 X 坐标（北向坐标，单位：米）
     * <p>
     * 在高斯-克吕格投影坐标系中，X 坐标表示点位在南北方向上的距离（北向坐标），
     * 数值上近似等于该点到赤道的子午线弧长（加上投影带假北距，通常为 0）。
     * 该坐标由 WGS84 纬度经横轴墨卡托投影正算公式转换得到，
     * 用于平面几何中的纵向距离计算和几何图形构建。
     * </p>
     * <p>
     * <b>取值范围：</b>理论上为 0 至约 10,000,000 米（对应北纬 0° 至 90°），
     * 实际取值取决于投影带和点位纬度。
     * </p>
     * <p>
     * <b>生成方式：</b>通过 GeoTools 的 {@link org.opengis.referencing.operation.MathTransform}
     * 将 WGS84 坐标正向转换得到，转换过程涉及椭球参数、投影带中央经线和假东距等参数。
     * </p>
     */
    private double gaussX;
    /**
     * 高斯投影 Y 坐标（东向坐标，单位：米）
     * <p>
     * 在高斯-克吕格投影坐标系中，Y 坐标表示点位在东西方向上的距离（东向坐标），
     * 以投影带中央经线为基准，向东为正、向西为负，并叠加假东距（False Easting）
     * 以确保所有 Y 坐标均为正值（便于工程应用）。
     * 该坐标由 WGS84 经度经横轴墨卡托投影正算公式转换得到，
     * 用于平面几何中的横向距离计算和几何图形构建。
     * </p>
     * <p>
     * <b>取值特点：</b>中央经线处 Y 值等于假东距（通常为 500,000 米 + 带号 × 1,000,000 米），
     * 向东西两侧对称递减/递增。
     * </p>
     * <p>
     * <b>生成方式：</b>与 {@link #gaussX} 同时通过坐标转换生成，共享同一套投影参数。
     * </p>
     */
    private double gaussY;

    /**
     * 点位所属多边形索引
     * <p>
     * 用于标识该点在聚类或空间分析后归属于哪个多边形（地块）。
     * 在 DBSCAN 聚类、道路拆分或多边形合并等算法中，
     * 该字段记录点位的聚类簇编号或地块序号，便于后续按多边形分组统计和溯源。
     * </p>
     * <p>
     * <b>取值说明：</b>
     * <ul>
     *   <li>null：表示该点尚未被分配到任何多边形（初始状态或噪声点）；</li>
     *   <li>≥ 0 的整数：表示该点所属的多边形索引，与聚类结果或拆分结果中的序号对应。</li>
     * </ul>
     * </p>
     */
    private Integer polygonIndex;

    /**
     * 默认无参构造方法
     * <p>
     * 创建一个空的 GaussPoint 实例，所有坐标值为 0.0，polygonIndex 为 null。
     * 通常用于反射创建、序列化反序列化或后续通过 setter 方法填充数据。
     * </p>
     */
    public GaussPoint() {
    }

    /**
     * 获取高斯投影 X 坐标（北向坐标，单位：米）
     *
     * @return 高斯投影 X 坐标（米），若未设置则返回 0.0
     */
    public double getGaussX() {
        return gaussX;
    }

    /**
     * 设置高斯投影 X 坐标（北向坐标，单位：米）
     * <p>
     * 通常由坐标转换工具方法自动设置，不建议手动修改，
     * 以确保与 WGS84 经纬度坐标和投影参数的一致性。
     * </p>
     *
     * @param gaussX 高斯投影 X 坐标（米），必须大于等于 0
     */
    public void setGaussX(double gaussX) {
        this.gaussX = gaussX;
    }

    /**
     * 获取高斯投影 Y 坐标（东向坐标，单位：米）
     *
     * @return 高斯投影 Y 坐标（米），若未设置则返回 0.0
     */
    public double getGaussY() {
        return gaussY;
    }

    /**
     * 设置高斯投影 Y 坐标（东向坐标，单位：米）
     * <p>
     * 通常由坐标转换工具方法自动设置，不建议手动修改，
     * 以确保与 WGS84 经纬度坐标和投影参数的一致性。
     * </p>
     *
     * @param gaussY 高斯投影 Y 坐标（米），通常为正值（含假东距）
     */
    public void setGaussY(double gaussY) {
        this.gaussY = gaussY;
    }

    /**
     * 获取点位所属多边形索引
     *
     * @return 多边形索引（≥ 0 的整数），若尚未分配则返回 null
     */
    public Integer getPolygonIndex() {
        return polygonIndex;
    }

    /**
     * 设置点位所属多边形索引
     * <p>
     * 通常在聚类或空间分析完成后由算法自动设置，
     * 用于标识该点在哪个多边形（地块）内。
     * </p>
     *
     * @param polygonIndex 多边形索引（≥ 0 的整数），null 表示未分配或噪声点
     */
    public void setPolygonIndex(Integer polygonIndex) {
        this.polygonIndex = polygonIndex;
    }
}
