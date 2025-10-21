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
    <version>1.0_jdk8_x64</version>
    <classifier>shaded</classifier>
</dependency>
```

### 使用示例，更多调用方法请看JtsUtil.java源码

```java

```

```java

```