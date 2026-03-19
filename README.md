# JTS工具类

## 功能特点

1. **坐标系统转换**：WGS84与高斯-克吕格投影坐标系双向转换，支持批量点转换和几何图形转换，具备线程安全的坐标转换器缓存机制

2. **空间几何分析**：
   - 点-几何关系判断（点是否在圆形区域内）
   - 多边形交集计算，支持WKT格式几何的相交轮廓提取
   - 几何缓冲区生成和合并操作

3. **距离与面积计算**：
   - 基于Haversine算法的球面距离计算
   - 球面面积计算（支持多边形和多多边形，含孔洞处理）
   - 国家标准亩数换算（1亩=666.6667平方米，保留4位小数精度）

4. **轨迹数据处理**：
   - 智能轨迹点过滤（零坐标、无效时间、地理范围验证）
   - 轨迹分段处理（基于距离阈值和时间间隔）
   - Douglas-Peucker抽稀优化
   - 分块并行处理支持万级轨迹点

5. **农业应用专用**：
   - 地块创建：基于轨迹点和作业宽度生成作业区域
   - 道路拆分：智能识别并分离轨迹中的道路部分
   - 作业面积精确计算，支持农业补贴核算

6. **高性能架构**：
   - 多级缓存机制（CRS缓存、坐标转换器缓存）
   - STRtree空间索引优化
   - 内存预分配和智能简化策略
   - 线程安全设计，支持高并发场景

## 环境

* jdk8 x64 及以上版本

## 依赖

```xml

<dependency>
    <groupId>sunyu.util</groupId>
    <artifactId>util-gis</artifactId>
   <!-- {util.version}_{jdk.version} -->
   <version>3.5_jdk8</version>
    <classifier>shaded</classifier>
</dependency>
```


## 测试工具

### 聚类分析调试工具

项目提供两种聚类测试工具，支持出差现场调试：

1. **Java桌面版** (`ClusterVisualizationGui.java`)：基于Swing的图形界面工具，支持高斯投影坐标文件加载、DBSCAN聚类参数配置、可视化结果展示

2. **Python网页版** (`ClusterVisualizationGui.py`)：基于Dash框架的Web可视化工具，支持文件上传、交互式参数调节、实时聚类结果展示

### 功能对比

| 功能项 | Java桌面版 | Python网页版 |
| ------ | ---------- | ------------ |
| 界面交互 | 本地窗口 | Web浏览器 |
| 数据加载 | 本地文件 | 本地文件 |
| 聚类参数 | 手动输入 | 滑动条调节 |
| 实时展示 | 本地窗口 | 浏览器页面 |
| 系统要求 | Windows/Linux/Mac | 任何操作系统 |

### 几何轮廓与相交测试工具

1. **WKT相交测试** (`showIntersection.html`)：基于天地图API的Web测试工具，支持输入两个WKT几何轮廓，可视化展示相交重叠情况，使用红蓝两色区分不同几何体，支持多种几何类型（点、线、面、多面）

2. **综合几何测试** (`showGeometrysTemplate_leaflet.html`)：基于Leaflet和天地图的综合测试平台，功能包括：
   - 地块轮廓展示与亩数计算
   - 轨迹数据可视化与回放
   - 点位标注与信息展示
   - 距离测量与面积测量
   - 多几何体批量展示
   - 实时坐标显示与比例尺
   - 支持轨迹播放控制与速度调节

### 注意事项

- 聚类分析工具均支持eps和minPts参数实时调节，可视化展示聚类结果，便于现场分析轨迹数据的空间分布特征
- 几何测试工具基于Web地图API，支持本地文件环境使用，无需额外安装，直接在浏览器中打开即可使用
- 所有测试工具均提供直观的可视化界面，支持实时交互和参数调节，适合现场调试和数据分析


## 使用方式

#### 工具生命周期管理
```java
// 构建器模式创建工具实例
GisUtil gisUtil = GisUtil.builder().build();

// 释放资源（坐标转换器缓存、CRS缓存等）
gisUtil.close();
```

#### 轨迹数据处理
```java
// 过滤异常轨迹点（零坐标、无效时间、地理范围验证）
List<Wgs84Point> filteredPoints = gisUtil.filterWgs84Points(wgs84Points);
```

#### 距离与面积计算
```java
// Haversine算法计算两点间球面距离（米）
double distance = gisUtil.haversine(wgs84Point1, wgs84Point2);

// 计算WGS84几何图形的球面面积（平方米）
double area = gisUtil.calculateSphericalArea(wgs84Geometry);

// 计算WGS84几何图形的面积（亩，保留4位小数）
double mu = gisUtil.calcMu(wgs84Geometry);

// 根据WKT字符串计算面积（亩）
double mu = gisUtil.calcMu(wgs84Wkt);
```

#### 空间几何分析

```java
import sunyu.util.pojo.FarmPlot;

// 判断点是否在圆形区域内（中心点+半径）
boolean inCircle = gisUtil.inCircle(wgs84Point, wgs84CenterPoint, radius);

// 判断WGS84点是否在几何图形内
boolean inGeometry = gisUtil.inGeometry(wgs84Point, wgs84Geometry);

// 判断点是否在矩形区域内（左上角+右下角）
boolean inRectangle = gisUtil.inRectangle(wgs84Point, wgs84TopLeftPoint, wgs84BottomRightPoint);

// 计算两个WKT几何图形的相交轮廓
WktIntersectionResult result = gisUtil.intersection(wgs84WKT1, wgs84WKT2);
String intersectionWkt = result.getWkt();
double intersectionMu = result.getMu();

// 查找距离目标点最近的点（带容差）
Wgs84Point closestPoint = gisUtil.findClosestPoint(targetWgs84Point, wgs84Points, tolerance);

// 查找距离目标点列表最近的多个点
List<Wgs84Point> closestPoints = gisUtil.findClosestPointList(targetPointList, wgs84Points);

// 合并多个WKT
String mergeWkt = gisUtil.mergeWgs84WKTStr(wgs84WktList);

// 合并多个WKT为一个FarmPlot对象，里面包含亩数
FarmPlot mergeWkt = gisUtil.mergeWgs84WKT(wgs84WktList);

// 将WKT字符串转换为4维数组
double[][][][] wktTo4DArray = gisUtil.wktTo4DArray(wgs84Wkt);
```

#### 坐标系统转换
```java
// WGS84点列表转高斯投影点列表
List<GaussPoint> gaussPoints = gisUtil.toGaussPointList(wgs84Points);

// 高斯投影点列表转WGS84点列表  
List<Wgs84Point> wgs84Points = gisUtil.toWgs84PointList(gaussPoints);

// WGS84几何转高斯投影几何
Geometry gaussGeometry = gisUtil.toGaussGeometry(wgs84Geometry);

// 高斯投影几何转WGS84几何
Geometry wgs84Geometry = gisUtil.toWgs84Geometry(gaussGeometry);

// WKT字符串转WGS84几何
Geometry geometry = gisUtil.toWgs84Geometry(wgs84WKT);
```

#### 农业应用专用
```java
// 根据轨迹点和作业宽度创建地块信息
FarmPlot farmPlot = gisUtil.getFarmPlot(wgs84Points, workingWidth);

// 智能作业轨迹道路拆分（默认参数）
SplitResult splitResult = gisUtil.splitRoad(wgs84Points, workingWidth);

// 智能作业轨迹道路拆分（自定义参数）
SplitResult splitResult = gisUtil.splitRoad(wgs84Points, workingWidth, splitRoadParams);
```