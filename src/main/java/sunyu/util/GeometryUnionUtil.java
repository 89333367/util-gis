package sunyu.util;

import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.GeometryCollection;
import org.locationtech.jts.geom.GeometryFactory;

import java.util.ArrayList;
import java.util.List;

/**
 * 几何图形合并工具类
 * 提供高效的几何图形合并方法
 *
 * @author SunYu
 */
public class GeometryUnionUtil {

    private final GeometryFactory geometryFactory;

    public GeometryUnionUtil(GeometryFactory geometryFactory) {
        this.geometryFactory = geometryFactory;
    }

    /**
     * 高效合并几何图形列表
     * 使用分组批量合并策略避免大量单独union操作的性能问题
     *
     * @param geometries 待合并的几何图形列表
     *
     * @return 合并后的几何图形
     */
    public Geometry union(List<Geometry> geometries) {
        if (geometries == null || geometries.isEmpty()) {
            return null;
        }

        if (geometries.size() == 1) {
            return geometries.get(0);
        }

        // 对于小量数据，直接使用GeometryCollection.union()
        if (geometries.size() <= 1000) {
            GeometryCollection collection = geometryFactory.createGeometryCollection(
                    geometries.toArray(new Geometry[0]));
            return collection.union();
        }

        // 对于大量数据，使用分组合并策略
        return unionInGroups(geometries, 100); // 每组100个
    }

    /**
     * 分组合并几何图形
     *
     * @param geometries 待合并的几何图形列表
     * @param groupSize  每组的大小
     *
     * @return 合并后的几何图形
     */
    private Geometry unionInGroups(List<Geometry> geometries, int groupSize) {
        if (geometries.size() <= groupSize) {
            GeometryCollection collection = geometryFactory.createGeometryCollection(
                    geometries.toArray(new Geometry[0]));
            return collection.union();
        }

        // 分组处理
        List<Geometry> groups = new ArrayList<>();
        for (int i = 0; i < geometries.size(); i += groupSize) {
            int end = Math.min(i + groupSize, geometries.size());
            List<Geometry> group = geometries.subList(i, end);

            GeometryCollection collection = geometryFactory.createGeometryCollection(
                    group.toArray(new Geometry[0]));
            groups.add(collection.union());
        }

        // 递归合并组
        return unionInGroups(groups, groupSize);
    }
}