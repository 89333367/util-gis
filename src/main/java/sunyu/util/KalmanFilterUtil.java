package sunyu.util;

import org.geotools.geometry.jts.JTS;
import org.geotools.referencing.CRS;
import org.geotools.referencing.crs.DefaultGeographicCRS;
import org.locationtech.jts.geom.Coordinate;
import org.opengis.referencing.crs.CoordinateReferenceSystem;
import org.opengis.referencing.operation.MathTransform;
import sunyu.util.pojo.TrackPoint;

import java.time.Duration;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.stream.Collectors;

/**
 * 卡尔曼滤波轨迹处理工具类
 * 核心功能：对GPS等轨迹数据进行滤波去噪，自动计算速度/方向角（不依赖外部传入），修正异常点（如GPS跳变）
 * 核心原理：基于匀速运动模型的卡尔曼滤波，通过"预测-更新"循环平衡轨迹平滑性与响应性
 * 使用注意：
 * 1. 非线程安全，多线程环境需创建独立实例调用filterTrack方法
 * 2. 依赖GeoTools库实现坐标转换（经纬度↔墨卡托），需确保依赖正确加载
 * 3. 适用于车辆、行人等中低速运动轨迹，默认参数已优化常规场景
 *
 * @author SunYu
 */
public class KalmanFilterUtil {
    /**
     * WGS84坐标系（EPSG:4326）：全球通用的经纬度坐标系
     * 坐标单位：经度（°）、纬度（°），无法直接计算两点间直线距离，需转墨卡托投影
     */
    private static final CoordinateReferenceSystem WGS84_CRS = DefaultGeographicCRS.WGS84;

    /**
     * 墨卡托坐标系（EPSG:3857）：Web地图常用投影坐标系
     * 坐标单位：米（x：东向，y：北向），支持平面欧几里得距离计算，是滤波的核心坐标系
     */
    private static final CoordinateReferenceSystem MERCATOR_CRS;

    /**
     * 经纬度→墨卡托坐标转换工具：将WGS84经纬度转换为墨卡托平面坐标，用于距离和速度计算
     */
    private static MathTransform wgs84ToMercator;

    /**
     * 墨卡托→经纬度坐标转换工具：将滤波后的墨卡托坐标转回WGS84经纬度，用于最终输出
     */
    private static MathTransform mercatorToWgs84;

    /**
     * 卡尔曼滤波状态矩阵：存储当前时刻的核心状态
     * 维度：4行1列，元素含义：
     * [0][0]：x坐标（墨卡托，米，东向）
     * [1][0]：y坐标（墨卡托，米，北向）
     * [2][0]：vx速度分量（米/秒，东向，正为东，负为西）
     * [3][0]：vy速度分量（米/秒，北向，正为北，负为南）
     */
    private double[][] stateMatrix;

    /**
     * 误差协方差矩阵P：描述状态矩阵中每个元素的不确定性（误差大小）
     * 维度：4行4列，对角线元素分别对应x、y、vx、vy的误差，非对角线为交叉误差
     * 初始值越大，代表初始状态的不确定性越高
     */
    private double[][] errorCovMatrix;

    /**
     * 测量噪声矩阵R：描述测量值（经纬度转墨卡托后的坐标）的噪声大小
     * 维度：2行2列（仅关注x、y坐标的测量噪声，速度不直接测量）
     * 对角线元素为测量噪声标准差的平方，值越大，代表越不信任测量值
     */
    private final double[][] measurementNoise;

    /**
     * 过程噪声标准差：描述运动过程中的不确定性（如加速度变化、突发转向）
     * 作用：值越大，滤波越信任测量值（响应快，平滑差）；值越小，越信任预测值（平滑好，响应慢）
     * 默认值5.0：平衡中低速运动（如车辆）的机动需求（如转弯）和平滑需求（去噪）
     */
    private double processNoiseStd = 5.0;

    /**
     * 测量噪声标准差：描述测量设备（如GPS）的误差大小
     * 作用：值越大，滤波越不信任测量值（减少异常点影响）；值越小，越信任测量值（易受噪声干扰）
     * 默认值50.0：针对GPS常见的10-50米误差，减少异常跳变对速度的影响
     */
    private double measurementNoiseStd = 50.0;

    /**
     * 残差限制阈值（米）：测量值与预测值的偏差（残差）超过此值时，强制裁剪到阈值
     * 作用：避免极端异常的测量偏差（如GPS跳变200米）导致速度计算异常（如瞬间几百km/h）
     * 默认值50米：参考GPS典型误差范围，兼顾去噪和轨迹真实性
     */
    private static final double RESIDUAL_LIMIT = 50.0;

    /**
     * 静态初始化块：初始化坐标系转换工具，仅执行一次
     * 逻辑：加载EPSG:3857墨卡托坐标系，创建经纬度与墨卡托的双向转换工具
     * 异常处理：若初始化失败（如GeoTools依赖缺失），抛出运行时异常，提示检查依赖
     */
    static {
        try {
            // 加载墨卡托坐标系（EPSG:3857），true表示强制使用指定CRS，不自动适配
            MERCATOR_CRS = CRS.decode("EPSG:3857", true);
            // 创建经纬度→墨卡托的转换工具（false表示不允许坐标系自动转换，确保精度）
            wgs84ToMercator = CRS.findMathTransform(WGS84_CRS, MERCATOR_CRS, false);
            // 创建墨卡托→经纬度的转换工具，用于最终输出
            mercatorToWgs84 = CRS.findMathTransform(MERCATOR_CRS, WGS84_CRS, false);
        } catch (Exception e) {
            // 坐标系初始化失败会导致后续所有转换异常，因此直接抛出运行时异常终止程序
            throw new RuntimeException("坐标系初始化失败，请检查GeoTools依赖（建议版本33.0+）", e);
        }
    }

    /**
     * 自定义噪声参数构造函数：适用于需要根据具体场景调整滤波特性的场景
     *
     * @param processNoiseStd     过程噪声标准差（推荐范围：3.0-10.0，机动大的场景（如无人机）设大值）
     * @param measurementNoiseStd 测量噪声标准差（推荐范围：30.0-100.0，GPS噪声大设大值）
     *                            逻辑：1. 赋值自定义噪声参数 2. 初始化测量噪声矩阵R（对角线为标准差的平方）3. 初始化误差协方差矩阵P
     */
    public KalmanFilterUtil(double processNoiseStd, double measurementNoiseStd) {
        // 赋值自定义过程噪声标准差，覆盖默认值
        this.processNoiseStd = processNoiseStd;
        // 赋值自定义测量噪声标准差，覆盖默认值
        this.measurementNoiseStd = measurementNoiseStd;
        // 初始化测量噪声矩阵R：2行2列对角矩阵，x和y方向噪声独立（值为标准差的平方）
        this.measurementNoise = new double[2][2];
        this.measurementNoise[0][0] = Math.pow(this.measurementNoiseStd, 2); // x方向测量噪声
        this.measurementNoise[1][1] = Math.pow(this.measurementNoiseStd, 2); // y方向测量噪声
        // 初始化误差协方差矩阵P，设置初始状态的不确定性
        initializeErrorCovMatrix();
    }

    /**
     * 默认构造函数：适用于大多数中低速轨迹场景（如车辆、行人）
     * 逻辑：调用自定义参数构造函数，传入默认噪声参数（processNoiseStd=5.0，measurementNoiseStd=50.0）
     * 适用场景：无需精细调参，追求平衡的去噪效果和轨迹响应速度
     */
    public KalmanFilterUtil() {
        this(5.0, 50.0); // 传入默认噪声参数，简化调用
    }

    /**
     * 初始化误差协方差矩阵P：设置初始状态的不确定性
     * 逻辑：对角线元素分别对应x、y、vx、vy的初始误差，非对角线为0（初始假设各状态独立）
     * 数值说明：
     * - x/y初始误差10米：假设初始位置有10米左右的不确定性（如GPS初始定位误差）
     * - vx/vy初始误差8米/秒：假设初始速度有8m/s（≈28.8km/h）的不确定性，避免初始速度异常
     */
    private void initializeErrorCovMatrix() {
        // 初始化4行4列的误差协方差矩阵
        errorCovMatrix = new double[4][4];
        errorCovMatrix[0][0] = 10.0;  // x坐标初始误差（米）
        errorCovMatrix[1][1] = 10.0;  // y坐标初始误差（米）
        errorCovMatrix[2][2] = 8.0;   // vx速度初始误差（米/秒）
        errorCovMatrix[3][3] = 8.0;   // vy速度初始误差（米/秒）
    }

    /**
     * 对外暴露的预处理方法：获取计算好速度/方向角的轨迹点列表
     * 作用：供外部验证预处理结果（如查看速度是否正确计算），避免直接依赖内部方法
     *
     * @param originalPoints 原始轨迹点列表（仅需包含经纬度和时间，速度/方向角可选）
     *
     * @return 预处理后的轨迹点列表（速度/方向角已自动计算，基于墨卡托距离）
     */
    public List<TrackPoint> getPreprocessedTrack(List<TrackPoint> originalPoints) {
        // 调用内部预处理方法，返回结果
        return preprocessTrackPoints(originalPoints);
    }

    /**
     * 内部轨迹预处理方法：核心作用是清洗原始数据，计算速度/方向角，为滤波做准备
     * 预处理步骤：1. 过滤无效点 2. 按时间排序 3. 转换墨卡托坐标 4. 计算速度/方向角
     * 关键设计：生成新的TrackPoint对象，不修改原始对象，避免引用混乱导致的数据污染
     *
     * @param originalPoints 原始轨迹点列表（可能包含无效点、乱序、无速度/方向角）
     *
     * @return 预处理后的轨迹点列表（有效、有序、含正确速度/方向角）
     */
    private List<TrackPoint> preprocessTrackPoints(List<TrackPoint> originalPoints) {
        // 1. 过滤无效点+按时间排序+生成新对象：
        // - 过滤time为null的点（无时间无法计算速度和时序）
        // - 按时间升序排序（确保轨迹是按运动时序排列，避免逆序导致速度为负）
        // - map生成新TrackPoint对象（不修改原始对象，避免外部调用者数据被篡改）
        List<TrackPoint> sortedNewPoints = originalPoints.stream()
                .filter(point -> point.getTime() != null) // 过滤无时间的无效点
                .sorted(Comparator.comparing(TrackPoint::getTime)) // 按时间升序排序
                .map(point -> new TrackPoint(point.getLon(), point.getLat(), point.getTime())) // 生成新对象
                .collect(Collectors.toList());

        // 若预处理后不足2个点（无法计算速度和滤波），直接返回（后续滤波会跳过）
        if (sortedNewPoints.size() < 2) {
            return sortedNewPoints;
        }

        // 2. 预计算所有点的墨卡托坐标：
        // - 原因：经纬度无法直接计算直线距离（地球是球体），墨卡托投影后可按平面距离计算
        // - 预计算：避免后续循环中重复转换，提升效率
        List<Coordinate> mercatorCoords = new ArrayList<>();
        for (TrackPoint point : sortedNewPoints) {
            // 将当前点的经纬度转换为墨卡托坐标
            Coordinate mercator = convertToMercator(point.getLon(), point.getLat());
            // 加入墨卡托坐标列表，与轨迹点一一对应
            mercatorCoords.add(mercator);
        }

        // 3. 遍历计算每个点的速度和方向角：
        // - 第1个点：无前置点，用与第2个点的关系反向计算（确保所有点都有速度/方向角）
        // - 非第1个点：用与前一个点的关系计算（基于时序的连续运动）
        for (int i = 0; i < sortedNewPoints.size(); i++) {
            // 当前处理的轨迹点（新对象）
            TrackPoint curr = sortedNewPoints.get(i);
            // 当前点的墨卡托坐标（用于距离计算）
            Coordinate currMercator = mercatorCoords.get(i);
            // 初始化速度（km/h）和方向角（°）
            double speedKmh = 0.0;
            double direction = 0.0;

            // 处理第1个点：无前置点，用第2个点计算
            if (i == 0) {
                // 获取第2个轨迹点（作为当前点的"下一个点"）
                TrackPoint next = sortedNewPoints.get(1);
                // 获取第2个点的墨卡托坐标
                Coordinate nextMercator = mercatorCoords.get(1);
                // 计算当前点与第2个点的时间差（秒）：Duration.between获取毫秒差，转秒
                double dt = Duration.between(curr.getTime(), next.getTime()).toMillis() / 1000.0;
                // 计算两点间墨卡托距离（米）：欧几里得距离公式√[(x2-x1)²+(y2-y1)²]
                double distanceM = Math.sqrt(Math.pow(nextMercator.x - currMercator.x, 2) + Math.pow(nextMercator.y - currMercator.y, 2));
                // 计算速度（km/h）：速度=距离/时间（m/s）→ 乘以3.6转换为km/h（日常速度单位）
                speedKmh = (distanceM / dt) * 3.6;
                // 计算方向角（°）：基于墨卡托坐标的东向和北向偏移，正北为0°顺时针递增
                direction = calculateDirection(currMercator, nextMercator);
                // 打印日志，调试确认速度计算是否正确
            } else {
                // 处理非第1个点：用前一个点计算
                // 获取前一个轨迹点
                TrackPoint prev = sortedNewPoints.get(i - 1);
                // 获取前一个点的墨卡托坐标
                Coordinate prevMercator = mercatorCoords.get(i - 1);
                // 计算当前点与前一个点的时间差（秒）
                double dt = Duration.between(prev.getTime(), curr.getTime()).toMillis() / 1000.0;
                // 计算两点间墨卡托距离（米）：欧几里得公式
                double distanceM = Math.sqrt(Math.pow(currMercator.x - prevMercator.x, 2) + Math.pow(currMercator.y - prevMercator.y, 2));
                // 计算速度（km/h）：距离/时间 → 转km/h
                speedKmh = (distanceM / dt) * 3.6;
                // 计算方向角（°）：基于前一个点到当前点的偏移
                direction = calculateDirection(prevMercator, currMercator);
                // 打印日志，调试确认速度计算是否正确
            }

            // 将计算好的速度和方向角设置到当前新轨迹点（覆盖默认0值）
            curr.setSpeed(speedKmh);
            curr.setDirection(direction);
        }

        // 返回预处理后的轨迹点列表（有效、有序、含速度/方向角）
        return sortedNewPoints;
    }

    /**
     * 轨迹滤波核心方法：输入原始轨迹点，输出滤波去噪后的轨迹点
     * 核心逻辑：基于匀速运动模型的卡尔曼滤波，通过"预测→更新"循环修正轨迹，处理异常点
     * 关键特性：不依赖外部传入的速度/方向角，完全基于经纬度和时间自动计算
     *
     * @param points 原始轨迹点列表（需包含经纬度和时间，速度/方向角可选）
     *
     * @return 滤波后的轨迹点列表（含去噪后的位置、速度、方向角）
     *
     * @throws IllegalArgumentException 若输入列表为null或包含无效点（time为null）
     */
    public List<TrackPoint> filterTrack(List<TrackPoint> points) {
        // 输入校验1：防止输入null列表导致后续stream操作空指针
        if (points == null) {
            throw new IllegalArgumentException("输入轨迹点列表不能为null，请传入有效轨迹数据");
        }

        // 输入校验2：防止列表中包含null点或time为null的点（导致时间计算异常）
        for (TrackPoint point : points) {
            if (point == null) {
                throw new IllegalArgumentException("轨迹点列表中不能包含null元素，请检查原始数据");
            }
            if (point.getTime() == null) {
                throw new IllegalArgumentException("轨迹点必须设置时间（time不能为null），请补充时间信息");
            }
        }

        // 步骤1：对原始轨迹进行预处理（过滤、排序、计算速度/方向角）
        List<TrackPoint> preprocessed = preprocessTrackPoints(points);
        // 若预处理后不足2个点（无法进行滤波，因为需要至少2个点计算速度），直接返回
        if (preprocessed.size() < 2) {
            return new ArrayList<>(preprocessed);
        }

        // 步骤2：初始化滤波结果列表，存储最终滤波后的轨迹点
        List<TrackPoint> filtered = new ArrayList<>();
        // 获取预处理后的第1个点（作为滤波的初始状态）
        TrackPoint first = preprocessed.get(0);
        // 将第1个点的经纬度转换为墨卡托坐标（用于状态矩阵初始化）
        Coordinate firstMercator = convertToMercator(first.getLon(), first.getLat());
        // 将第1个点的速度从km/h转换为m/s（状态矩阵中速度单位为m/s）
        double firstSpeedMs = kmhToMs(first.getSpeed());
        // 将第1个点的方向角从度转换为弧度（三角函数计算需弧度）
        double firstDirRadian = Math.toRadians(first.getDirection());

        // 步骤3：初始化卡尔曼滤波的状态矩阵（初始状态基于第1个预处理点）
        stateMatrix = new double[4][1];
        stateMatrix[0][0] = firstMercator.x; // 初始x坐标（墨卡托，米）
        stateMatrix[1][0] = firstMercator.y; // 初始y坐标（墨卡托，米）
        // 初始vx速度分量（东向）：速度×sin(方向角)，sin(θ)对应x轴（东向）投影
        stateMatrix[2][0] = firstSpeedMs * Math.sin(firstDirRadian);
        // 初始vy速度分量（北向）：速度×cos(方向角)，cos(θ)对应y轴（北向）投影
        stateMatrix[3][0] = firstSpeedMs * Math.cos(firstDirRadian);

        // 步骤4：将第1个预处理点加入滤波结果（第1个点无前置滤波过程，直接用预处理值）
        filtered.add(new TrackPoint(
                first.getLon(), // 原始经纬度（无需滤波，因是初始点）
                first.getLat(),
                first.getTime(),
                first.getSpeed(), // 预处理计算的速度
                first.getDirection() // 预处理计算的方向角
        ));

        // 步骤5：遍历处理后续所有预处理点（从第2个点开始，索引i=1）
        for (int i = 1; i < preprocessed.size(); i++) {
            // 当前处理的预处理点（含正确的速度/方向角，但位置可能有异常）
            TrackPoint curr = preprocessed.get(i);
            // 前一个预处理点（用于计算时间差）
            TrackPoint prev = preprocessed.get(i - 1);
            // 计算当前点与前一个点的时间差（秒）：避免时间差过小导致速度异常
            double dt = Duration.between(prev.getTime(), curr.getTime()).toMillis() / 1000.0;
            // 限制时间差≥0.1秒：防止时间差趋近于0导致速度趋近于无穷大（除以0异常）
            if (dt <= 0.1) dt = 0.1;

            // 卡尔曼滤波第一步：预测（基于前一状态预测当前状态）
            predict(dt);

            // 将当前预处理点的经纬度转换为墨卡托坐标（作为测量值，用于更新）
            Coordinate currMercator = convertToMercator(curr.getLon(), curr.getLat());
            // 卡尔曼滤波第二步：更新（用测量值修正预测状态，带残差限制避免异常）
            updateWithResidualLimit(currMercator.x, currMercator.y);

            // 步骤6：生成当前点的滤波结果（墨卡托坐标转经纬度，计算滤波后速度/方向角）
            // 将滤波后的墨卡托坐标转回WGS84经纬度（用户需要的最终坐标格式）
            Coordinate filteredWgs84 = convertToWgs84(stateMatrix[0][0], stateMatrix[1][0]);
            // 计算滤波后的速度（m/s）：合成vx和vy的合速度，√(vx²+vy²)
            double filteredSpeedMs = Math.sqrt(Math.pow(stateMatrix[2][0], 2) + Math.pow(stateMatrix[3][0], 2));
            // 限制速度上限为初始速度的2倍：防止极端异常（如滤波后速度仍过大），兼顾合理性
            filteredSpeedMs = Math.min(filteredSpeedMs, firstSpeedMs * 2);
            // 计算滤波后的方向角（°）：基于vx和vy的方向，atan2(vx, vy)→弧度转度，调整到0~360°
            double filteredDirDegree = (Math.toDegrees(Math.atan2(stateMatrix[2][0], stateMatrix[3][0])) + 360) % 360;

            // 将滤波结果加入列表：包含去噪后的经纬度、时间、速度、方向角
            filtered.add(new TrackPoint(
                    filteredWgs84.x, // 滤波后的经度
                    filteredWgs84.y, // 滤波后的纬度
                    curr.getTime(),  // 保留原始时间（时间无需滤波）
                    msToKmh(filteredSpeedMs), // 速度从m/s转换为km/h
                    filteredDirDegree // 滤波后的方向角（0~360°）
            ));
        }

        // 返回最终的滤波轨迹点列表
        return filtered;
    }

    /**
     * 带残差限制的卡尔曼更新步骤：核心是用测量值修正预测状态，同时限制异常偏差
     * 与普通更新的区别：对测量值与预测值的偏差（残差）进行裁剪，避免异常偏差导致速度突变
     *
     * @param measuredX 测量的x坐标（当前点的墨卡托坐标，米）
     * @param measuredY 测量的y坐标（当前点的墨卡托坐标，米）
     */
    private void updateWithResidualLimit(double measuredX, double measuredY) {
        // 1. 定义测量矩阵H：描述"状态→测量值"的映射关系
        // 维度：2行4列，仅提取状态矩阵中的x（[0][0]）和y（[1][0]）作为测量值（速度不直接测量）
        double[][] measurementMatrix = {{1, 0, 0, 0}, {0, 1, 0, 0}};
        // 计算测量矩阵H的转置矩阵H_T：矩阵运算需要（卡尔曼增益公式要求）
        double[][] H_T = transposeMatrix(measurementMatrix);

        // 2. 计算卡尔曼增益K：权衡预测值和测量值的权重
        // 公式步骤：
        // S = H * P * H_T + R （测量预测误差=预测误差投影+测量噪声）
        // K = P * H_T * S⁻¹ （增益=预测误差投影 / 测量预测误差，值越大越信任测量）
        double[][] S_temp1 = multiplyMatrix(measurementMatrix, errorCovMatrix); // H*P
        double[][] S_temp2 = multiplyMatrix(S_temp1, H_T); // H*P*H_T
        double[][] S = addMatrix(S_temp2, measurementNoise); // S = H*P*H_T + R
        double[][] K_temp = multiplyMatrix(errorCovMatrix, H_T); // P*H_T
        double[][] K = multiplyMatrix(K_temp, invert2x2Matrix(S)); // K = P*H_T * S⁻¹（S⁻¹是S的逆矩阵）

        // 3. 计算残差：测量值与预测测量值的偏差（反映测量与预测的差异）
        // z：测量值向量（2行1列，x和y）
        double[][] z = {{measuredX}, {measuredY}};
        // H_x：预测的测量值（基于前一状态预测的当前x和y）
        double[][] H_x = multiplyMatrix(measurementMatrix, stateMatrix);
        // residual = z - H_x （偏差=实际测量 - 预测测量）
        double[][] residual = subtractMatrix(z, H_x);

        // 4. 残差裁剪：限制偏差在±RESIDUAL_LIMIT（50米）内，避免极端异常
        // 原因：若测量值跳变200米，残差200米会导致状态更新过大，速度瞬间几百km/h
        residual[0][0] = Math.max(Math.min(residual[0][0], RESIDUAL_LIMIT), -RESIDUAL_LIMIT); // x方向残差裁剪
        residual[1][0] = Math.max(Math.min(residual[1][0], RESIDUAL_LIMIT), -RESIDUAL_LIMIT); // y方向残差裁剪

        // 5. 更新状态矩阵：用裁剪后的残差修正预测状态
        // x_new = x_pred + K * residual （新状态=预测状态 + 增益*偏差，平衡预测和测量）
        stateMatrix = addMatrix(stateMatrix, multiplyMatrix(K, residual));

        // 6. 更新误差协方差矩阵：反映新状态的不确定性
        // P_new = (I - K*H) * P_pred （新误差=预测误差 * 增益修正系数，确保误差非负）
        double[][] I = {{1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}}; // 4阶单位矩阵
        double[][] K_H = multiplyMatrix(K, measurementMatrix); // K*H（增益与测量矩阵的乘积）
        errorCovMatrix = multiplyMatrix(subtractMatrix(I, K_H), errorCovMatrix); // P_new = (I-KH)*P_pred
    }

    /**
     * 卡尔曼预测步骤：基于前一时刻的状态，预测当前时刻的状态
     * 核心假设：匀速运动模型（无加速度），适用于中低速运动（如车辆、行人）
     *
     * @param dt 前一时刻到当前时刻的时间差（秒）
     */
    private void predict(double dt) {
        // 1. 定义状态转移矩阵F：描述"前一状态→当前预测状态"的映射（匀速模型）
        // 维度：4行4列，公式依据：
        // x_new = x_old + vx_old * dt （当前x=前x + 前x方向速度*时间）
        // y_new = y_old + vy_old * dt （当前y=前y + 前y方向速度*时间）
        // vx_new = vx_old （匀速假设，x方向速度不变）
        // vy_new = vy_old （匀速假设，y方向速度不变）
        double[][] transitionMatrix = new double[4][4];
        transitionMatrix[0][0] = 1; // x_new = 1*x_old + 0*y_old + dt*vx_old + 0*vy_old
        transitionMatrix[1][1] = 1; // y_new = 0*x_old + 1*y_old + 0*vx_old + dt*vy_old
        transitionMatrix[2][2] = 1; // vx_new = 0*x_old + 0*y_old + 1*vx_old + 0*vy_old
        transitionMatrix[3][3] = 1; // vy_new = 0*x_old + 0*y_old + 0*vx_old + 1*vy_old
        transitionMatrix[0][2] = dt; // x与vx的关联系数（dt）
        transitionMatrix[1][3] = dt; // y与vy的关联系数（dt）

        // 2. 预测状态矩阵：基于转移矩阵计算当前预测状态
        // x_pred = F * x_old （当前预测状态=转移矩阵*前一状态）
        stateMatrix = multiplyMatrix(transitionMatrix, stateMatrix);

        // 3. 预测误差协方差矩阵：考虑过程噪声（运动不确定性）
        // 公式：P_pred = F * P_old * F_T + Q （预测误差=前误差转移 + 过程噪声）
        double[][] F_T = transposeMatrix(transitionMatrix); // F的转置矩阵
        double[][] P_temp = multiplyMatrix(transitionMatrix, errorCovMatrix); // F*P_old
        double[][] P_pred = multiplyMatrix(P_temp, F_T); // F*P_old*F_T
        // 获取过程噪声矩阵Q（基于时间差dt动态计算，匀速模型的噪声分布）
        double[][] Q = getProcessNoiseMatrix(dt);
        // 预测误差 = 转移后的误差 + 过程噪声
        errorCovMatrix = addMatrix(P_pred, Q);
    }

    /**
     * 计算过程噪声矩阵Q：描述运动过程中的不确定性（如加速度变化）
     * 矩阵维度：4行4列，基于匀速运动模型的噪声公式，与时间差dt相关
     *
     * @param dt 时间差（秒）：dt越大，过程噪声越大（运动不确定性越高）
     *
     * @return 4行4列的过程噪声矩阵Q
     */
    private double[][] getProcessNoiseMatrix(double dt) {
        // 计算dt的各次幂：用于过程噪声公式
        double dt2 = dt * dt; // dt²
        double dt3 = dt2 * dt; // dt³
        double dt4 = dt3 * dt; // dt⁴
        // 过程噪声功率：过程噪声标准差的平方（噪声大小的量化）
        double q = Math.pow(processNoiseStd, 2);

        // 初始化过程噪声矩阵Q：基于匀速模型的噪声分布
        // 公式依据：匀速运动中，位置误差与dt⁴成正比，速度误差与dt²成正比
        double[][] processNoise = new double[4][4];
        processNoise[0][0] = dt4 / 4 * q; // x位置的过程噪声
        processNoise[0][2] = dt3 / 2 * q; // x与vx的交叉噪声
        processNoise[1][1] = dt4 / 4 * q; // y位置的过程噪声
        processNoise[1][3] = dt3 / 2 * q; // y与vy的交叉噪声
        processNoise[2][0] = dt3 / 2 * q; // vx与x的交叉噪声
        processNoise[2][2] = dt2 * q;    // vx速度的过程噪声
        processNoise[3][1] = dt3 / 2 * q; // vy与y的交叉噪声
        processNoise[3][3] = dt2 * q;    // vy速度的过程噪声

        return processNoise;
    }

    /**
     * 计算两点间的方向角：基于墨卡托坐标，正北为0°，顺时针递增（0°北、90°东、180°南、270°西）
     *
     * @param prevMercator 前一个点的墨卡托坐标（起点）
     * @param currMercator 当前点的墨卡托坐标（终点）
     *
     * @return 方向角（°），范围0~360°
     */
    private double calculateDirection(Coordinate prevMercator, Coordinate currMercator) {
        // 计算东向偏移量（x轴差）：curr.x - prev.x（正为东，负为西）
        double deltaX = currMercator.x - prevMercator.x;
        // 计算北向偏移量（y轴差）：curr.y - prev.y（正为北，负为南）
        double deltaY = currMercator.y - prevMercator.y;
        // 计算方向角弧度：atan2(deltaX, deltaY) → 结果范围-π~π（对应-180°~180°）
        // 原理：atan2(对边, 邻边)，此处对边为东向偏移，邻边为北向偏移，对应正北顺时针角度
        double radian = Math.atan2(deltaX, deltaY);
        // 弧度转角度：-180°~180°
        double degree = Math.toDegrees(radian);
        // 调整角度范围到0~360°：若为负角，加360°（如-90°→270°）
        return (degree + 360) % 360;
    }

    /**
     * WGS84经纬度转换为墨卡托坐标（EPSG:3857）
     *
     * @param lon 经度（°）：范围-180~180，东经为正，西经为负
     * @param lat 纬度（°）：范围-85.0511~85.0511（墨卡托投影的有效纬度范围）
     *
     * @return 墨卡托坐标（x：东向米，y：北向米）
     *
     * @throws RuntimeException 若转换失败（如纬度超出有效范围）
     */
    private Coordinate convertToMercator(double lon, double lat) {
        try {
            // 创建WGS84经纬度坐标对象（x=经度，y=纬度）
            Coordinate wgs84 = new Coordinate(lon, lat);
            // 创建墨卡托坐标对象，用于存储转换结果
            Coordinate mercator = new Coordinate();
            // 执行转换：将经纬度坐标转换为墨卡托坐标
            JTS.transform(wgs84, mercator, wgs84ToMercator);
            // 返回转换后的墨卡托坐标
            return mercator;
        } catch (Exception e) {
            // 转换失败（如纬度超出±85.0511°），抛出异常提示具体经纬度
            throw new RuntimeException("经纬度(" + lon + "," + lat + ")转墨卡托失败，可能纬度超出有效范围（±85.0511°）", e);
        }
    }

    /**
     * 墨卡托坐标（EPSG:3857）转换为WGS84经纬度
     *
     * @param x 墨卡托x坐标（米，东向）
     * @param y 墨卡托y坐标（米，北向）
     *
     * @return WGS84经纬度坐标（x=经度°，y=纬度°）
     *
     * @throws RuntimeException 若转换失败（如坐标超出墨卡托投影范围）
     */
    private Coordinate convertToWgs84(double x, double y) {
        try {
            // 创建墨卡托坐标对象（x=东向米，y=北向米）
            Coordinate mercator = new Coordinate(x, y);
            // 创建WGS84经纬度坐标对象，用于存储转换结果
            Coordinate wgs84 = new Coordinate();
            // 执行转换：将墨卡托坐标转换为经纬度
            JTS.transform(mercator, wgs84, mercatorToWgs84);
            // 返回转换后的经纬度坐标
            return wgs84;
        } catch (Exception e) {
            // 转换失败（如墨卡托坐标超出合理范围），抛出异常提示具体坐标
            throw new RuntimeException("墨卡托(" + x + "," + y + ")转经纬度失败，可能坐标超出投影范围", e);
        }
    }

    /**
     * 速度单位转换：千米/小时（km/h）→ 米/秒（m/s）
     * 公式：1 km/h = 1000 m / 3600 s = 1/3.6 m/s
     *
     * @param kmh 速度（km/h）
     *
     * @return 速度（m/s）
     */
    private double kmhToMs(double kmh) {
        return kmh / 3.6;
    }

    /**
     * 速度单位转换：米/秒（m/s）→ 千米/小时（km/h）
     * 公式：1 m/s = 3.6 km/h（如10 m/s = 36 km/h）
     *
     * @param ms 速度（m/s）
     *
     * @return 速度（km/h）
     */
    private double msToKmh(double ms) {
        return ms * 3.6;
    }

    /**
     * 矩阵乘法：计算矩阵A × 矩阵B，要求A的列数 = B的行数
     *
     * @param a 矩阵A（rowsA行 × colsA列）
     * @param b 矩阵B（colsA行 × colsB列）
     *
     * @return 结果矩阵C（rowsA行 × colsB列）
     */
    private double[][] multiplyMatrix(double[][] a, double[][] b) {
        // 获取矩阵A的行数
        int rowsA = a.length;
        // 获取矩阵A的列数（假设A非空，取第一行的长度）
        int colsA = a[0].length;
        // 获取矩阵B的列数（假设B非空，取第一行的长度）
        int colsB = b[0].length;

        // 初始化结果矩阵C（rowsA行 × colsB列），初始值为0
        double[][] result = new double[rowsA][colsB];

        // 三重循环计算矩阵乘积：C[i][j] = Σ(A[i][k] × B[k][j])（k从0到colsA-1）
        for (int i = 0; i < rowsA; i++) {
            for (int k = 0; k < colsA; k++) {
                // 优化：若A[i][k]为0，跳过当前k（减少乘法运算，提升效率）
                if (a[i][k] == 0) continue;
                for (int j = 0; j < colsB; j++) {
                    result[i][j] += a[i][k] * b[k][j];
                }
            }
        }

        return result;
    }

    /**
     * 矩阵转置：将矩阵rows行×cols列 转换为 cols行×rows列
     *
     * @param matrix 原始矩阵（rows行 × cols列）
     *
     * @return 转置后的矩阵（cols行 × rows列）
     */
    private double[][] transposeMatrix(double[][] matrix) {
        // 获取原始矩阵的行数
        int rows = matrix.length;
        // 获取原始矩阵的列数（假设矩阵非空）
        int cols = matrix[0].length;

        // 初始化转置矩阵（cols行 × rows列）
        double[][] transposed = new double[cols][rows];

        // 循环赋值：转置矩阵的[j][i] = 原始矩阵的[i][j]
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                transposed[j][i] = matrix[i][j];
            }
        }

        return transposed;
    }

    /**
     * 矩阵加法：计算矩阵A + 矩阵B，要求A和B的行数、列数完全相同
     *
     * @param a 矩阵A（rows行 × cols列）
     * @param b 矩阵B（rows行 × cols列）
     *
     * @return 结果矩阵C（rows行 × cols列），C[i][j] = A[i][j] + B[i][j]
     */
    private double[][] addMatrix(double[][] a, double[][] b) {
        // 获取矩阵A的行数
        int rows = a.length;
        // 获取矩阵A的列数（假设A非空）
        int cols = a[0].length;

        // 初始化结果矩阵C（rows行 × cols列）
        double[][] result = new double[rows][cols];

        // 循环计算：对应元素相加
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                result[i][j] = a[i][j] + b[i][j];
            }
        }

        return result;
    }

    /**
     * 矩阵减法：计算矩阵A - 矩阵B，要求A和B的行数、列数完全相同
     *
     * @param a 矩阵A（被减矩阵，rows行 × cols列）
     * @param b 矩阵B（减矩阵，rows行 × cols列）
     *
     * @return 结果矩阵C（rows行 × cols列），C[i][j] = A[i][j] - B[i][j]
     */
    private double[][] subtractMatrix(double[][] a, double[][] b) {
        // 获取矩阵A的行数
        int rows = a.length;
        // 获取矩阵A的列数（假设A非空）
        int cols = a[0].length;

        // 初始化结果矩阵C（rows行 × cols列）
        double[][] result = new double[rows][cols];

        // 循环计算：对应元素相减
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                result[i][j] = a[i][j] - b[i][j];
            }
        }

        return result;
    }

    /**
     * 2×2矩阵求逆：计算2阶矩阵的逆矩阵，仅用于卡尔曼滤波中的S矩阵（2×2）求逆
     * 逆矩阵存在条件：行列式det ≠ 0（det = a*d - b*c，矩阵为[[a,b],[c,d]]）
     *
     * @param matrix 2×2输入矩阵（[[a,b],[c,d]]）
     *
     * @return 2×2逆矩阵（[[d/det, -b/det], [-c/det, a/det]]）
     *
     * @throws ArithmeticException 若矩阵不可逆（det接近0）
     */
    private double[][] invert2x2Matrix(double[][] matrix) {
        // 计算矩阵的行列式det = a*d - b*c（matrix[0][0]=a, matrix[0][1]=b, matrix[1][0]=c, matrix[1][1]=d）
        double det = matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];

        // 检查行列式是否接近0（小于1e-8）：若接近0，矩阵不可逆，抛出异常
        if (Math.abs(det) < 1e-8) {
            throw new ArithmeticException("2×2矩阵不可逆，行列式=" + det + "，可能测量噪声设置过小");
        }

        // 计算行列式的倒数（用于逆矩阵公式）
        double invDet = 1.0 / det;

        // 初始化逆矩阵：基于2×2矩阵逆公式
        // 逆矩阵 = (1/det) × [[d, -b], [-c, a]]
        double[][] inverse = new double[2][2];
        inverse[0][0] = matrix[1][1] * invDet; // d/det
        inverse[0][1] = -matrix[0][1] * invDet; // -b/det
        inverse[1][0] = -matrix[1][0] * invDet; // -c/det
        inverse[1][1] = matrix[0][0] * invDet; // a/det

        return inverse;
    }
}