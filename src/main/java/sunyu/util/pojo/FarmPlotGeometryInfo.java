package sunyu.util.pojo;

import org.locationtech.jts.geom.Geometry;

import java.util.List;
import java.util.Map;

public class FarmPlotGeometryInfo {
    private Map<Integer, Geometry> geometryMap;
    private Map<Integer, List<GaussPoint>> geometryGaussPointMap;

    public FarmPlotGeometryInfo(Map<Integer, Geometry> geometryMap, Map<Integer, List<GaussPoint>> geometryGaussPointMap) {
        this.geometryMap = geometryMap;
        this.geometryGaussPointMap = geometryGaussPointMap;
    }


    public Map<Integer, Geometry> getGeometryMap() {
        return geometryMap;
    }

    public void setGeometryMap(Map<Integer, Geometry> geometryMap) {
        this.geometryMap = geometryMap;
    }

    public Map<Integer, List<GaussPoint>> getGeometryGaussPointMap() {
        return geometryGaussPointMap;
    }

    public void setGeometryGaussPointMap(Map<Integer, List<GaussPoint>> geometryGaussPointMap) {
        this.geometryGaussPointMap = geometryGaussPointMap;
    }
}
