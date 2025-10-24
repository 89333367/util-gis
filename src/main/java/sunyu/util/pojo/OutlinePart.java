package sunyu.util.pojo;

import java.time.LocalDateTime;

import org.locationtech.jts.geom.Geometry;

/**
 * 轮廓区块详情：单个 Polygon 的元数据。
 * @author SunYu
 */
public class OutlinePart {
    /**
     * 区块的几何图形（Polygon）
     */
    private Geometry polygon;
    /**
     * 区块的起始时间
     */
    private LocalDateTime startTime;
    /**
     * 区块的结束时间
     */
    private LocalDateTime endTime;
    /**
     * 区块面积（亩）
     */
    private double mu;
    /**
     * 区块的WKT表示
     */
    private String wkt;

    public OutlinePart(Geometry polygon, LocalDateTime startTime, LocalDateTime endTime, double mu, String wkt) {
        this.polygon = polygon;
        this.startTime = startTime;
        this.endTime = endTime;
        this.mu = mu;
        this.wkt = wkt;
    }

    public Geometry getPolygon() {
        return polygon;
    }

    public LocalDateTime getStartTime() {
        return startTime;
    }

    public LocalDateTime getEndTime() {
        return endTime;
    }

    public double getMu() {
        // 为确保与 Turf.js 一致，这里按 WGS84 球面面积计算并转换为亩
        this.mu = calcMuSpherical(this.polygon);
        return mu;
    }

    // 使用 Turf.js 的球面环面积公式计算 Polygon 的亩数（保留三位小数）
    private static double calcMuSpherical(org.locationtech.jts.geom.Geometry geom) {
        if (!(geom instanceof org.locationtech.jts.geom.Polygon)) {
            // 非 Polygon，返回构造时的 mu 或 0
            return 0.0;
        }
        org.locationtech.jts.geom.Polygon p = (org.locationtech.jts.geom.Polygon) geom;
        double areaOuter = ringArea(p.getExteriorRing());
        double holesArea = 0.0;
        for (int i = 0; i < p.getNumInteriorRing(); i++) {
            holesArea += ringArea(p.getInteriorRingN(i));
        }
        double areaSqm = Math.abs(areaOuter) - Math.abs(holesArea);
        // 1 亩 = 666.6667 平方米
        return Math.round((areaSqm / 666.6667) * 1000.0) / 1000.0;
    }

    // 球面环面积（平方米），与 Turf.js 的 ringArea 保持一致
    private static double ringArea(org.locationtech.jts.geom.LineString ring) {
        org.locationtech.jts.geom.Coordinate[] coords = ring.getCoordinates();
        int len = (coords == null) ? 0 : coords.length;
        if (len <= 2)
            return 0.0;
        double area = 0.0;
        for (int i = 0; i < len; i++) {
            org.locationtech.jts.geom.Coordinate p1, p2, p3;
            if (i == len - 2) {
                p1 = coords[i];
                p2 = coords[i + 1];
                p3 = coords[0];
            } else if (i == len - 1) {
                p1 = coords[i];
                p2 = coords[0];
                p3 = coords[1];
            } else {
                p1 = coords[i];
                p2 = coords[i + 1];
                p3 = coords[i + 2];
            }
            area += (toRad(p3.x) - toRad(p1.x)) * Math.sin(toRad(p2.y));
        }
        double R = 6378137.0; // Turf 默认采用的半径（WGS84）
        return area * R * R / 2.0;
    }

    private static double toRad(double deg) {
        return deg * Math.PI / 180.0;
    }

    public String getWkt() {
        return wkt;
    }
}
