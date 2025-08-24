# JTS工具类

## 功能特点

1. 拥有电子围栏判断功能
2. 离散点生成几何图形功能
3. 集合图形计算面积功能
4. geojson.html用于在浏览器中展示WKT信息
5. map.html用于在地图上画点、线、矩形、圆形、多边形，F12查看经纬度信息

## 还可以判断如下几何图形关系

- 相等(Equals)：几何形状拓扑上相等。
- 脱节(Disjoint)：几何形状没有共有的点。
- 相交(Intersects)：几何形状至少有一个共有点（区别于脱节）
- 接触(Touches)：几何形状有至少一个公共的边界点，但是没有内部点。
- 交叉(Crosses)：几何形状共享一些但不是所有的内部点。
- 内含(Within)：几何形状A的线都在几何形状B内部。
- 包含(Contains)：几何形状B的线都在几何形状A内部（区别于内含）
- 重叠(Overlaps)：几何形状共享一部分但不是所有的公共点，而且相交处有他们自己相同的区域。

## 环境

* jdk8 x64 及以上版本

## 依赖

```xml

<dependency>
    <groupId>sunyu.util</groupId>
    <artifactId>util-gis</artifactId>
    <!-- {util.version}_{jdk.version}_{architecture.version} -->
    <version>1.0_jdk11_x64</version>
    <classifier>shaded</classifier>
</dependency>
```

### 使用示例，更多调用方法请看JtsUtil.java源码

```java
public class TestJtsUtil {
    Log log = LogFactory.get();
    JtsUtil jtsUtil = JtsUtil.builder().build();

    @Test
    public void 测试多边形面积() throws ParseException {
        String wkt = ResourceUtil.readUtf8Str("polygon.txt");
        log.debug(wkt);
        double area = jtsUtil.geometryArea(jtsUtil.getWktReader().read(wkt));
        log.debug("多边形面积：{} 平方米", area);
        //平方米换为亩，计算口诀为“加半左移三”
        //1平方米＝0.0015亩，如128平方米等于多少亩?计算方法是先用128加128的一半：128＋64＝192，再把小数点左移3位，即得出亩数为0.192
        log.debug("多边形面积：{} 亩", NumberUtil.div(area + NumberUtil.div(area, 2), 1000));
        log.debug("多边形面积：{} 亩", area / 666.67);
    }

    @Test
    public void 测试点是否在矩形内() {
        boolean out = jtsUtil.inRectangle(116.30699, 40.05496, 116.291323, 40.057832, 116.303828, 40.051536);
        Assert.isFalse(out);
        log.debug("是否在矩形内：{}", out);
        boolean in = jtsUtil.inRectangle(116.301456, 40.055623, 116.291323, 40.057832, 116.303828, 40.051536);
        Assert.isTrue(in);
        log.debug("是否在矩形内：{}", in);
    }


    @Test
    public void 测试点是否在圆内() {
        boolean out = jtsUtil.inCircle(116.311158, 40.055043, 116.302247, 40.05797, 427.25225281516754);
        Assert.isFalse(out);
        log.debug("是否在圆内：{}", out);
        boolean in = jtsUtil.inCircle(116.304474, 40.055319, 116.302247, 40.05797, 427.25225281516754);
        Assert.isTrue(in);
        log.debug("是否在圆内：{}", in);
    }

    @Test
    public void 测试点是否在多边形内() throws ParseException {
        Geometry polygon = jtsUtil.getWktReader().read("POLYGON((116.312496 40.05944,116.314333 40.059426,116.314284 40.058046,116.312465 40.058129,116.313826 40.05876,116.312496 40.05944))");
        Geometry outPoint = jtsUtil.getWktReader().read("POINT(116.313484 40.058736)");
        Assert.isFalse(jtsUtil.inGeometry((Point) outPoint, polygon));
        log.debug("是否在多边形内：{}", jtsUtil.inGeometry((Point) outPoint, polygon));
        Geometry inPoint = jtsUtil.getWktReader().read("POINT(116.31299 40.058187)");
        Assert.isTrue(jtsUtil.inGeometry((Point) inPoint, polygon));
        log.debug("是否在多边形内：{}", jtsUtil.inGeometry((Point) inPoint, polygon));
    }
}
```

```java
public class TestCoordTransformUtil {
    Log log = LogFactory.get();

    @Test
    void t001() {
        CoordTransformUtil coordTransformUtil = CoordTransformUtil.builder().build();
        log.debug("{}", coordTransformUtil.outOfChina(39.178762, 118.546783));
    }
}
```