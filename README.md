# JTS工具类

## 功能特点

1. 拓扑关系判断：相交、相等、脱节、接触、交叉、内含、包含、重叠（WGS84，支持 `POLYGON/MULTIPOLYGON`）
2. 电子围栏：点是否在多边形内（含边界），支持快速裁剪与拓扑修复
3. WKT 解析与转换：统一 WGS84 输出，支持坐标识别与高斯-克吕格分带投影转换
4. 面积与度量：球面公式面积（亩数换算）、两点 Haversine 距离（米）
5. 相交计算：两 WKT 求交并返回相交几何与亩数
6. 轨迹轮廓：根据轨迹点和总宽度生成轮廓（Polygon），含形态学处理
7. 道路分段：按总宽度对轨迹分段，支持最大段数控制
8. 示例页面：`src/test/resources/showGeometry.html` 在浏览器中展示几何图形关系示例

## 还可以判断如下几何图形关系（DE-9IM 拓扑语义）

- 相等（Equals）：两几何拓扑结构相等（`equalsTopo`）
- 脱节（Disjoint）：没有任何公共点（`disjoint`）
- 相交（Intersects）：至少有一个公共点（`intersects`）
- 接触（Touches）：存在公共边界点但没有内部公共点（`touches`）
- 交叉（Crosses）：共享部分内部点，交集既非包含也非相等（`crosses`）
- 内含（Within）：A 完全在 B 内（`within`）
- 包含（Contains）：B 完全在 A 内（`contains`）
- 重叠（Overlaps）：同维几何部分重合，交集既不是包含也不是相等（`overlaps`）

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

### 使用示例，更多调用方法请看 `GisUtil.java` 源码

```java
// 拓扑关系与电子围栏
GisUtil gis = GisUtil.builder().build();

String wktA = "POLYGON((120 30, 120.01 30, 120.01 30.01, 120 30.01, 120 30))";
String wktB = "POLYGON((120.005 30.005, 120.015 30.005, 120.015 30.015, 120.005 30.015, 120.005 30.005))";

boolean intersects = gis.intersects(wktA, wktB);
boolean equalsTopo = gis.equalsWkt(wktA, wktB);
boolean disjoint = gis.disjoint(wktA, wktB);
boolean touches = gis.touches(wktA, wktB);
boolean crosses = gis.crosses(wktA, wktB);
boolean within = gis.within(wktA, wktB);
boolean contains = gis.contains(wktA, wktB);
boolean overlaps = gis.overlaps(wktA, wktB);

boolean inFence = gis.pointInPolygon(new CoordinatePoint(120.008, 30.008), wktA);
System.out.printf("intersects=%s, inFence=%s\n", intersects, inFence);
```

```java
// 相交计算与面积、几何转换
GisUtil gis = GisUtil.builder().build();

WktIntersectionResult res = gis.intersection(wktA, wktB);
System.out.println("inter WKT=" + res.getWkt());
System.out.println("inter mu=" + res.getMu());

Geometry g = gis.fromWkt(wktA); // 解析并转换到高斯-克吕格米制坐标
String wkt = gis.toWkt(g);      // 统一输出为 WGS84 WKT

// 两点大圆距离（米）
double d = gis.haversine(new CoordinatePoint(120.0, 30.0), new CoordinatePoint(120.01, 30.02));
System.out.println("distance(m)=" + d);
```

## 轨迹轮廓与分段示例

```java
import java.time.LocalDateTime;
import java.util.Arrays;
import java.util.List;
import sunyu.util.GisUtil;
import sunyu.util.pojo.TrackPoint;
import sunyu.util.pojo.OutlinePart;
import sunyu.util.pojo.SplitRoadResult;

// 准备轨迹点（WGS84，经纬度，时间）
List<TrackPoint> seg = Arrays.asList(
    new TrackPoint(LocalDateTime.parse("2024-05-01T08:00:00"), 120.0000, 30.0000),
    new TrackPoint(LocalDateTime.parse("2024-05-01T08:05:00"), 120.0010, 30.0005),
    new TrackPoint(LocalDateTime.parse("2024-05-01T08:10:00"), 120.0020, 30.0010),
    new TrackPoint(LocalDateTime.parse("2024-05-01T08:15:00"), 120.0030, 30.0015)
);

GisUtil gis = GisUtil.builder().build();

// 1) 单体轮廓（不拆分）：总宽度（米，左右合计）
OutlinePart outline = gis.getOutline(seg, 12.0);
System.out.println("mu=" + outline.getMu());
System.out.println("wkt=" + outline.getWkt());
System.out.println("start=" + outline.getStartTime() + ", end=" + outline.getEndTime());
System.out.println("trackStr=" + outline.getTrackStr()); // "lon,lat,yyyyMMddHHmmss#..."

// 2) 分段轮廓（保留面积 Top-N）
SplitRoadResult split = gis.splitRoad(seg, 12.0, 5); // 最多保留5段
System.out.println("outlineWkt=" + split.getWkt());
for (OutlinePart part : split.getParts()) {
    System.out.println("part mu=" + part.getMu() + ", start=" + part.getStartTime() + ", end=" + part.getEndTime());
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

- `toWkt(Geometry)`：几何转 WKT（统一 WGS84，坐标合法性识别与修复）。
- `fromWkt(String)`：WKT 解析为 `Geometry` 并转换到高斯-克吕格米制坐标（按首点经度分带）。
- `haversine(CoordinatePoint, CoordinatePoint)`：两点大圆距离，单位米（WGS84）。
- `calcMu(Geometry)` / `calcMu(String)`：几何面积换算为亩（球面公式，半径 6378137）。
- `intersection(String, String)`：两段 WKT 求交，返回 `WktIntersectionResult`（包含 `wkt` 与 `mu`）。
- 拓扑判断（WKT，支持 `POLYGON/MULTIPOLYGON`）：
  - `intersects`、`equalsWkt`、`disjoint`、`touches`、`crosses`、`within`、`contains`、`overlaps`。
- 轨迹轮廓：
  - `getOutline(List<TrackPoint>, double totalWidthM)`：线缓冲构建单体轮廓并返回 `OutlinePart`。
  - `splitRoad(List<TrackPoint>, double totalWidthM[, Integer maxSegments])`：保留面积 Top-N 的分段结果，返回 `SplitRoadResult`。

## 坐标系与单位

- 输入/输出 WKT 均使用 `WGS84 (EPSG:4326)`；形态学与面积计算在高斯-克吕格米制投影下进行（6 度分带，按首点经度自动选择）。
- `totalWidthM` 为作业总宽度（米），内部按单侧宽度（`widthM = totalWidthM / 2`）进行线缓冲。
- 面积单位为“亩”（基于球面面积公式，半径 6378137，与 Turf.js 对齐）。

## 行为与默认值

- 轨迹预处理：
  - 排除越界与 `(0,0)` 点并按时间升序；`splitRoad` 会额外按速度过滤（默认范围 `[1, 20] km/h`）。
  - 单体轮廓 `getOutline` 不做道路拆分、不平滑；如生成 `MultiPolygon`，仅做轻微开运算保留缝隙并选面积最大面。
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
- 轨迹点最少数量：`getOutline/splitRoad` 最少需要 3 个有效点。
- 结果为 `MultiPolygon`：`splitRoad` 根据面积裁剪保留不超过 `maxSegments` 的区块；`getOutline` 返回单体面。