package sunyu.util;

import java.util.ArrayList;
import java.util.List;

import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.GeometryFactory;

import cn.hutool.log.Log;
import cn.hutool.log.LogFactory;
import sunyu.util.pojo.SplitRoadResult;
import sunyu.util.pojo.TrackPoint;

/**
 * GIS工具类，封装轨迹轮廓构建、道路分段、坐标转换、拓扑判断与面积计算。
 * 
 * 
 * @author SunYu
 */
public class GisUtil implements AutoCloseable {
    private final Log log = LogFactory.get();
    // 配置参数，包含各种常量和默认值
    private final Config config;

    /**
     * 创建构建器实例（Builder）。
     * 通过构建器配置参数后再构建 `GisUtil`，避免直接实例化带来的不完整/不一致配置。
     *
     * @return 用于配置并构建 `GisUtil` 的构建器
     */
    public static Builder builder() {
        return new Builder();
    }

    /**
     * 私有构造函数，通过 Builder 创建实例。
     * 保持不可直接实例化，集中初始化日志与配置引用。
     *
     * @param config 内部配置对象，包含常量、默认值与缓存实例
     */
    private GisUtil(Config config) {
        log.info("[构建{}] 开始", this.getClass().getSimpleName());
        log.info("[构建{}] 结束", this.getClass().getSimpleName());
        this.config = config;
    }

    /**
     * 内部配置类。
     * <p>
     * 持有常量、默认参数、线程安全缓存与几何工厂；供 {@link GisUtil} 使用。
     * 该类包含所有GIS处理过程中需要的配置参数，确保线程安全的参数访问。
     * </p>
     */
    private static class Config {
        // 几何工厂
        private final GeometryFactory geometryFactory = new GeometryFactory();
        // 空几何集合 - 使用空数组参数更明确和安全
        private final Geometry EMPTYGEOM = geometryFactory.createGeometryCollection(new Geometry[0]);
        // 空WKT字符串
        private final String EMPTY_WKT = EMPTYGEOM.toText();
    }

    /**
     * 构建器类，用于构建GisUtil实例。
     * <p>
     * 实现构建器模式，允许灵活配置和创建GisUtil对象，提供流式API进行参数设置。
     * </p>
     */
    public static class Builder {
        // 配置对象，包含各种常量和默认值
        private final Config config = new Config();

        /**
         * 构建GisUtil实例。
         * <p>
         * 使用预配置的参数创建GisUtil对象，完成所有必要的初始化工作。
         * </p>
         * 
         * @return 已初始化完成的GisUtil对象，可直接使用
         */
        public GisUtil build() {
            return new GisUtil(config);
        }

    }

    /**
     * 关闭资源（AutoCloseable）
     * 当前仅清理内部缓存；预留扩展以释放外部资源。
     */
    @Override
    public void close() {
        log.info("[销毁{}] 开始", this.getClass().getSimpleName());
        log.info("[销毁{}] 结束", this.getClass().getSimpleName());
    }

    public SplitRoadResult splitRoad(List<TrackPoint> wgs84Points, double totalWidthM) {
        SplitRoadResult result = new SplitRoadResult();
        result.setOutline(config.EMPTYGEOM);
        result.setWkt(config.EMPTY_WKT);

        log.debug("过滤时间为空、经纬度异常、经纬度为0、速度为0、方向角异常的轨迹点");
        List<TrackPoint> filteredWgs84Points = new ArrayList<>();
        for (TrackPoint point : wgs84Points) {
            if (point.getTime() != null// 定位时间不能为空
                    && point.getLon() >= -180 && point.getLon() <= 180// 经度必须在[-180,180]之间
                    && point.getLat() >= -90 && point.getLat() <= 90// 纬度必须在[-90,90]之间
                    && point.getLon() != 0 && point.getLat() != 0// 经度和纬度不能为0
                    && point.getSpeed() > 0// 速度必须大于0
                    && point.getDirection() >= 0 && point.getDirection() <= 360) {// 方向角必须在[0,360]之间
                filteredWgs84Points.add(point);
            } else {
                log.warn("被过滤点位信息，时间：{} 经度：{} 纬度：{} 速度：{} 方向角：{}", point.getTime(), point.getLon(), point.getLat(),
                        point.getSpeed(), point.getDirection());
            }
        }

        // 按时间正序排序
        filteredWgs84Points.sort((p1, p2) -> p1.getTime().compareTo(p2.getTime()));

        return result;
    }

}