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