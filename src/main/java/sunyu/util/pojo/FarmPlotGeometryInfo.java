package sunyu.util.pojo;

import org.locationtech.jts.geom.Geometry;

import java.util.List;
import java.util.Map;

/**
 * 地块几何信息聚合实体类
 * <p>
 * 用于封装聚类或拆分后产生的多个地块几何图形的聚合信息。
 * 该类将每个地块的 JTS {@link Geometry} 对象与其对应的高斯投影轨迹点列表进行关联存储，
 * 以 {@link Integer} 类型的聚类索引（或地块序号）作为键，便于后续的空间分析、面积计算和轨迹溯源。
 * </p>
 * <p>
 * <b>数据结构说明：</b>
 * <ul>
 *   <li>{@link #geometryMap}：存储每个地块的高斯投影几何图形（JTS Geometry），键为聚类索引，值为几何对象；</li>
 *   <li>{@link #geometryGaussPointMap}：存储构成每个地块的原始高斯投影轨迹点集合，键与 geometryMap 对应，值为该地块内的所有轨迹点。</li>
 * </ul>
 * </p>
 * <p>
 * <b>典型使用场景：</b>
 * <ul>
 *   <li>DBSCAN 聚类后，将每个聚类簇的几何图形和原始点集合打包传递；</li>
 *   <li>道路拆分算法中，作为中间结果存储多个子地块的几何与轨迹信息；</li>
 *   <li>批量面积计算或 WKT 转换前的数据准备容器。</li>
 * </ul>
 * </p>
 *
 * @author SunYu
 * @see Geometry
 * @see GaussPoint
 * @since 1.0.0
 */
public class FarmPlotGeometryInfo {
    /**
     * 地块几何图形映射表
     * <p>
     * 键（Integer）：聚类索引或地块序号，从 0 开始递增；
     * 值（Geometry）：对应地块的高斯投影几何图形（Polygon 或 MultiPolygon），
     * 由轨迹点经缓冲、合并等空间操作生成，用于面积计算和 WKT 输出。
     * </p>
     */
    private Map<Integer, Geometry> geometryMap;
    /**
     * 地块轨迹点映射表
     * <p>
     * 键（Integer）：与 {@link #geometryMap} 中的键一一对应，确保几何图形与轨迹点的关联一致性；
     * 值（List&lt;GaussPoint&gt;）：构成该地块的所有原始高斯投影轨迹点，按时间顺序排列，
     * 支持作业里程计算、时间范围提取和轨迹回放等后续分析。
     * </p>
     */
    private Map<Integer, List<GaussPoint>> geometryGaussPointMap;

    /**
     * 全参构造方法
     * <p>
     * 同时初始化几何图形映射和轨迹点映射，要求两个映射的键集合保持一致，
     * 即每个聚类索引在两张表中均有对应数据。
     * </p>
     *
     * @param geometryMap           地块几何图形映射表（高斯投影坐标系），不可为 null
     * @param geometryGaussPointMap 地块轨迹点映射表（高斯投影坐标系），不可为 null
     */
    public FarmPlotGeometryInfo(Map<Integer, Geometry> geometryMap, Map<Integer, List<GaussPoint>> geometryGaussPointMap) {
        this.geometryMap = geometryMap;
        this.geometryGaussPointMap = geometryGaussPointMap;
    }


    /**
     * 获取地块几何图形映射表
     *
     * @return 地块几何图形映射表（高斯投影坐标系），键为聚类索引，值为 JTS Geometry 对象；
     * 若未初始化则可能返回 null
     */
    public Map<Integer, Geometry> getGeometryMap() {
        return geometryMap;
    }

    /**
     * 设置地块几何图形映射表
     * <p>
     * 通常与 {@link #setGeometryGaussPointMap(Map)} 同时调用，确保两个映射的键集合一致。
     * </p>
     *
     * @param geometryMap 地块几何图形映射表（高斯投影坐标系），键为聚类索引，值为 JTS Geometry 对象
     */
    public void setGeometryMap(Map<Integer, Geometry> geometryMap) {
        this.geometryMap = geometryMap;
    }

    /**
     * 获取地块轨迹点映射表
     *
     * @return 地块轨迹点映射表（高斯投影坐标系），键为聚类索引，值为按时间排序的 GaussPoint 列表；
     * 若未初始化则可能返回 null
     */
    public Map<Integer, List<GaussPoint>> getGeometryGaussPointMap() {
        return geometryGaussPointMap;
    }

    /**
     * 设置地块轨迹点映射表
     * <p>
     * 通常与 {@link #setGeometryMap(Map)} 同时调用，确保两个映射的键集合一致。
     * </p>
     *
     * @param geometryGaussPointMap 地块轨迹点映射表（高斯投影坐标系），键为聚类索引，值为 GaussPoint 列表
     */
    public void setGeometryGaussPointMap(Map<Integer, List<GaussPoint>> geometryGaussPointMap) {
        this.geometryGaussPointMap = geometryGaussPointMap;
    }
}
