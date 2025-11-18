# JTS工具类

## 功能特点

1. 坐标系统转换：WGS84与高斯-克吕格投影坐标系双向转换，支持批量点转换和几何图形转换
2. 空间几何分析：点-几何关系判断（点是否在圆内、矩形内、多边形内），支持多边形交集计算
3. 距离与面积计算：基于球面公式的精确距离计算（Haversine算法）和面积计算，支持亩数换算
4. 轨迹处理：智能轨迹分段处理，包含数据预处理、速度过滤、Douglas-Peucker抽稀优化、分块处理
5. 高性能优化：多级缓存机制、并行分块处理、智能空间索引优化，适合大数据量场景
6. 几何构建：基于轨迹点和作业宽度生成缓冲区轮廓，支持复杂轨迹的智能分段和合并

## 环境

* jdk8 x64 及以上版本

## 依赖

```xml

<dependency>
    <groupId>sunyu.util</groupId>
    <artifactId>util-gis</artifactId>
    <!-- {util.version}_{jdk.version}_{architecture.version} -->
    <version>1.0_jdk8_x64</version>
    <classifier>shaded</classifier>
</dependency>
```

### 使用示例

```java
import java.time.LocalDateTime;
import java.util.Arrays;
import java.util.List;
import sunyu.util.GisUtil;
import sunyu.util.pojo.TrackPoint;
import sunyu.util.pojo.CoordinatePoint;
import sunyu.util.pojo.SplitRoadResult;
import sunyu.util.pojo.OutlinePart;
import sunyu.util.pojo.WktIntersectionResult;
import org.locationtech.jts.geom.Geometry;

GisUtil gis = GisUtil.builder().build();

// 1) 坐标转换示例
TrackPoint wgs84Point = new TrackPoint(LocalDateTime.now(), 120.0, 30.0);
TrackPoint gaussPoint = gis.toGaussPoint(wgs84Point); // WGS84转高斯投影
System.out.println("WGS84: lon=" + wgs84Point.getLon() + ", lat=" + wgs84Point.getLat());
System.out.println("Gauss: x=" + gaussPoint.getLon() + ", y=" + gaussPoint.getLat());

// 2) 批量坐标转换
List<TrackPoint> wgs84Points = Arrays.asList(
    new TrackPoint(LocalDateTime.now(), 120.0, 30.0),
    new TrackPoint(LocalDateTime.now(), 120.01, 30.01)
);
List<TrackPoint> gaussPoints = gis.toGaussPointList(wgs84Points);
System.out.println("converted " + gaussPoints.size() + " points");

// 3) WKT字符串转几何对象
String wktPolygon = "POLYGON((120 30, 120.01 30, 120.01 30.01, 120 30.01, 120 30))";
Geometry geometry = gis.toWgs84Geometry(wktPolygon);
System.out.println("Geometry type: " + geometry.getGeometryType());

// 4) 距离计算（Haversine算法）
CoordinatePoint point1 = new CoordinatePoint(120.0, 30.0);
CoordinatePoint point2 = new CoordinatePoint(120.01, 30.02);
double distance = gis.haversine(point1, point2);
System.out.println("distance(m)=" + distance);

// 5) 点与几何关系判断
CoordinatePoint testPoint = new CoordinatePoint(120.005, 30.005);
boolean inPolygon = gis.inGeometry(testPoint, geometry);
System.out.println("point in polygon=" + inPolygon);

// 6) 圆形区域判断
CoordinatePoint center = new CoordinatePoint(120.0, 30.0);
boolean inCircle = gis.inCircle(testPoint, center, 1000.0); // 1公里半径
System.out.println("point in circle=" + inCircle);

// 7) 矩形区域判断
CoordinatePoint topLeft = new CoordinatePoint(119.99, 30.01);
CoordinatePoint bottomRight = new CoordinatePoint(120.01, 29.99);
boolean inRectangle = gis.inRectangle(testPoint, topLeft, bottomRight);
System.out.println("point in rectangle=" + inRectangle);
```

```java
// 面积计算与轨迹处理
GisUtil gis = GisUtil.builder().build();

// 1) 几何面积计算（亩数）
String wktA = "POLYGON((120 30, 120.01 30, 120.01 30.01, 120 30.01, 120 30))";
double mu = gis.calcMu(wktA);
System.out.println("area(mu)=" + mu);

// 2) 球面面积计算（平方米）
Geometry geometryA = gis.toWgs84Geometry(wktA);
double area = gis.calculateSphericalArea(geometryA);
System.out.println("area(sqm)=" + area);

// 3) 多边形交集计算
String wktB = "POLYGON((120.005 30.005, 120.015 30.005, 120.015 30.015, 120.005 30.015, 120.005 30.005))";
WktIntersectionResult intersection = gis.intersection(wktA, wktB);
System.out.println("intersection WKT=" + intersection.getWkt());
System.out.println("intersection mu=" + intersection.getMu());

// 4) 坐标系统转换 - 高斯转WGS84
Geometry gaussGeometry = gis.toGaussGeometry(geometryA);
Geometry backToWgs84 = gis.toWgs84Geometry(gaussGeometry);
System.out.println("converted back to WGS84: " + backToWgs84.getGeometryType());

// 5) 轨迹分段处理（农业作业轨迹）
List<TrackPoint> trackPoints = Arrays.asList(
    new TrackPoint(LocalDateTime.parse("2024-05-01T08:00:00"), 120.0000, 30.0000),
    new TrackPoint(LocalDateTime.parse("2024-05-01T08:05:00"), 120.0010, 30.0005),
    new TrackPoint(LocalDateTime.parse("2024-05-01T08:10:00"), 120.0020, 30.0010),
    new TrackPoint(LocalDateTime.parse("2024-05-01T08:15:00"), 120.0030, 30.0015)
);

// 智能轨迹分段，作业宽度12米
SplitRoadResult splitResult = gis.splitRoad(trackPoints, 12.0);
System.out.println("total area(mu)=" + splitResult.getMu());
System.out.println("total segments=" + splitResult.getParts().size());

for (OutlinePart part : splitResult.getParts()) {
    System.out.println("segment: " + part.getStartTime() + " ~ " + part.getEndTime() + 
                       ", area=" + String.format("%.2f", part.getMu()) + "mu");
    System.out.println("  WKT=" + part.getWkt());
}

// 6) 资源管理 - 释放GisUtil占用的资源
gis.close();
```

## 轨迹处理示例

```java
import java.time.LocalDateTime;
import java.util.Arrays;
import java.util.List;
import sunyu.util.GisUtil;
import sunyu.util.pojo.TrackPoint;
import sunyu.util.pojo.OutlinePart;
import sunyu.util.pojo.SplitRoadResult;

// 准备轨迹点（WGS84，经纬度，时间）
List<TrackPoint> trackPoints = Arrays.asList(
    new TrackPoint(LocalDateTime.parse("2024-05-01T08:00:00"), 120.0000, 30.0000),
    new TrackPoint(LocalDateTime.parse("2024-05-01T08:05:00"), 120.0010, 30.0005),
    new TrackPoint(LocalDateTime.parse("2024-05-01T08:10:00"), 120.0020, 30.0010),
    new TrackPoint(LocalDateTime.parse("2024-05-01T08:15:00"), 120.0030, 30.0015),
    new TrackPoint(LocalDateTime.parse("2024-05-01T08:20:00"), 120.0040, 30.0020)
);

GisUtil gis = GisUtil.builder().build();

// 智能轨迹分段处理：作业宽度12米
SplitRoadResult splitResult = gis.splitRoad(trackPoints, 12.0);
System.out.println("total area(mu)=" + splitResult.getMu());
System.out.println("total segments=" + splitResult.getParts().size());

for (OutlinePart part : splitResult.getParts()) {
    System.out.println("segment: " + part.getStartTime() + " ~ " + part.getEndTime() + 
                       ", area=" + String.format("%.2f", part.getMu()) + "mu");
    System.out.println("  WKT=" + part.getWkt());
}
```

## 数据结构

- `TrackPoint(time, lon, lat)`：轨迹点，时间为 `LocalDateTime`，经纬度为 WGS84。
- `OutlinePart`：单个区块的详情。
  - 关键字段：`getOutline()`、`getStartTime()`、`getEndTime()`、`getMu()`、`getWkt()`、`getTrackPoints()`、`getTrackStr()`、`getTotalWidthM()`。
- `SplitRoadResult`：分段结果汇总。
  - 关键字段：`getOutline()`、`getWkt()`、`getParts()`（按 `startTime` 升序返回）、`getTotalWidthM()`。
- `WktIntersectionResult`：相交结果，字段：`getWkt()`、`getMu()`。

## API速查

- 坐标系统转换：
  - `toGaussPoint(TrackPoint)`：单个WGS84坐标点转高斯投影坐标
  - `toGaussPointList(List<TrackPoint>)`：批量WGS84坐标点转高斯投影坐标
  - `toWgs84Geometry(Geometry)`：高斯投影几何图形转WGS84几何图形
  - `toWgs84Geometry(String)`：WGS84 WKT字符串转几何对象
  - `toGaussGeometry(Geometry)`：WGS84几何图形转高斯投影几何图形

- 距离与面积计算：
  - `haversine(CoordinatePoint, CoordinatePoint)`：两点大圆距离，单位米（WGS84）
  - `calcMu(Geometry)` / `calcMu(String)`：几何面积换算为亩（球面公式，半径 6378137）
  - `calculateSphericalArea(Geometry)`：计算几何图形的球面面积（平方米）

- 空间几何关系判断：
  - `inCircle(CoordinatePoint, CoordinatePoint, double)`：点是否在指定半径的圆形区域内
  - `inGeometry(CoordinatePoint, Geometry)`：点是否在多边形内（含边界）
  - `inRectangle(CoordinatePoint, CoordinatePoint, CoordinatePoint)`：点是否在矩形区域内

- 相交计算：
  - `intersection(String, String)`：两段WKT求交，返回 `WktIntersectionResult`（包含 `wkt` 与 `mu`）

- 轨迹处理：
  - `splitRoad(List<TrackPoint>, double totalWidthM)`：智能轨迹分段处理，返回 `SplitRoadResult`

## 坐标系与单位

- 输入/输出 WKT 均使用 `WGS84 (EPSG:4326)`；形态学与面积计算在高斯-克吕格米制投影下进行（6 度分带，按首点经度自动选择）。
- `totalWidthM` 为作业总宽度（米），内部按单侧宽度（`widthM = totalWidthM / 2`）进行线缓冲。
- 面积单位为“亩”（基于球面面积公式，半径 6378137，与 Turf.js 对齐）。

## 行为与默认值

- 轨迹预处理：
  - 排除越界与 `(0,0)` 点并按时间升序；`splitRoad` 会额外按速度过滤（默认范围 `[1, 20] km/h`）。
  - 单体轮廓处理不做道路拆分、不平滑；如生成 `MultiPolygon`，仅做轻微开运算保留缝隙并选面积最大面。
- 合法性修复：
  - 拓扑判断在发现非法坐标时会先执行 `buffer(0)` 修复（自相交等问题）。
- 可配置项：
  - 当前构建器未开放外部参数设置，内部默认值包含线简化公差、缓冲样式、缝隙轻微蚀刻、薄条裁剪等。

## 示例页面

- 打开 `src/test/resources/showGeometry.html`：在浏览器输入/粘贴 WKT（或轨迹），在天地图上查看轮廓与面积标注。
- 打开 `src/test/resources/showIntersection.html`：输入两段 WKT，查看彩色绘制与视图缩放效果。

## 构建与测试

- 构建库：
  - `mvn -q -DskipTests=true package`
- 运行测试（仅几何相关用例）：
  - `mvn -q -DskipTests=false test`
  - 注：少量数据库相关测试需要环境支持，若无可忽略这类用例。

## 常见问题

- WKT 自相交或坐标不合法：内部会尝试 `buffer(0)` 修复；建议输入合法几何。
- `POINT/LINESTRING` 支持：拓扑谓词以 `POLYGON/MULTIPOLYGON` 为主；其他类型用于展示与转换。
- 轨迹点最少数量：`splitRoad` 最少需要 3 个有效点。
- 结果为 `MultiPolygon`：`splitRoad` 根据面积裁剪保留不超过 `maxSegments` 的区块。