package sunyu.util;

import org.geotools.api.referencing.crs.CoordinateReferenceSystem;
import org.geotools.api.referencing.operation.MathTransform;
import org.geotools.geometry.jts.JTS;
import org.geotools.referencing.CRS;
import org.geotools.referencing.crs.DefaultGeographicCRS;
import org.locationtech.jts.geom.Coordinate;
import sunyu.util.pojo.TrackPoint;

import java.time.Duration;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.stream.Collectors;

/**
 * 卡尔曼滤波工具类
 */
public class KalmanFilterUtil {
    // 坐标系定义（x:东向，y:北向）
    private static final CoordinateReferenceSystem WGS84_CRS = DefaultGeographicCRS.WGS84;
    private static final CoordinateReferenceSystem MERCATOR_CRS;
    private static MathTransform wgs84ToMercator;
    private static MathTransform mercatorToWgs84;

    // 卡尔曼滤波核心参数
    private double[][] stateMatrix;       // [x(米,东), y(米,北), vx(米/秒,东向), vy(米/秒,北向)]
    private double[][] errorCovMatrix;    // 误差协方差矩阵 P
    private final double[][] measurementNoise;  // 测量噪声矩阵 R
    /**
     * 值较大时 (如 5.0-10.0)：
     * 滤波器更倾向于相信测量值而不是预测值
     * 对轨迹的快速变化响应更敏感
     * 平滑效果减弱，可能保留更多原始噪声
     * 适用于机动性强的目标跟踪
     * 值较小时 (如 0.1-1.0)：
     * 滤波器更相信预测值，平滑效果更强
     * 可能过度平滑，导致轨迹细节丢失
     * 对真实的位置变化响应较慢
     * 适用于平稳运动的目标
     */
    private double processNoiseStd = 2.0; // 过程噪声标准差
    /**
     * 值较大时 (如 10.0-20.0)：
     * 滤波器更相信预测值而非测量值
     * 平滑效果增强，但可能偏离真实轨迹
     * 对异常点和噪声的鲁棒性提高
     * 输出轨迹更平滑但可能滞后
     * 值较小时 (如 1.0-3.0)：
     * 滤波器更相信测量值
     * 对噪声敏感，可能导致轨迹抖动
     * 跟踪精度提高但稳定性下降
     * 快速响应真实位置变化
     */
    private double measurementNoiseStd = 8.0; // 测量噪声标准差

    static {
        try {
            MERCATOR_CRS = CRS.decode("EPSG:3857");
            wgs84ToMercator = CRS.findMathTransform(WGS84_CRS, MERCATOR_CRS, true);
            mercatorToWgs84 = CRS.findMathTransform(MERCATOR_CRS, WGS84_CRS, true);
        } catch (Exception e) {
            throw new RuntimeException("坐标系初始化失败", e);
        }
    }

    /**
     * // 高过程噪声 - 适用于频繁变向的车辆
     * new KalmanFilterUtil(5.0, 5.0);  // 对突然转向响应快
     * <p>
     * // 低过程噪声 - 适用于直线运动
     * new KalmanFilterUtil(0.5, 5.0);  // 强平滑效果，忽略小幅波动
     * <p>
     * // 高测量噪声 - 适用于低精度GPS
     * new KalmanFilterUtil(1.0, 15.0);  // 强噪声过滤，轨迹更平滑
     * <p>
     * // 低测量噪声 - 适用于高精度定位
     * new KalmanFilterUtil(1.0, 2.0);   // 紧跟原始轨迹，保留细节
     * <p>
     * // 默认平衡配置（适用于一般场景）
     * new KalmanFilterUtil(1.0, 5.0);
     * <p>
     * // 高平滑配置（适用于噪声大但运动平稳）
     * new KalmanFilterUtil(0.5, 10.0);
     * <p>
     * // 高响应配置（适用于机动性强的目标）
     * new KalmanFilterUtil(3.0, 3.0);
     * <p>
     * // 极度平滑配置（适用于高质量数据的精细处理）
     * new KalmanFilterUtil(0.1, 2.0);
     * <p>
     * // 城市道路车辆（频繁启停、变道）
     * new KalmanFilterUtil(3.0, 8.0);  // 高过程噪声，适应频繁机动
     * <p>
     * // 高速公路车辆（相对平稳运动）
     * new KalmanFilterUtil(1.5, 6.0);  // 中等过程噪声，平衡平滑与响应
     * <p>
     * // 公交车路线（固定路线，站点停靠）
     * new KalmanFilterUtil(2.0, 10.0); // 较高测量噪声，过滤站点停靠抖动
     * <p>
     * // 农田作业机械（直线作业，速度较慢）
     * new KalmanFilterUtil(0.8, 12.0); // 低过程噪声，高测量噪声
     * <p>
     * // 收割机作业（不规则路径，频繁转向）
     * new KalmanFilterUtil(2.5, 15.0); // 中等过程噪声，高测量噪声
     * <p>
     * // 播种机作业（相对规整的直线路径）
     * new KalmanFilterUtil(0.5, 8.0);  // 低过程噪声，中等测量噪声
     * <p>
     * // 无人机航拍（高度机动）
     * new KalmanFilterUtil(4.0, 6.0);  // 高过程噪声，适应快速机动
     * <p>
     * // 船舶航行（相对平稳的大范围移动）
     * new KalmanFilterUtil(1.0, 20.0); // 低过程噪声，高测量噪声
     * <p>
     * // 步行人员轨迹（频繁停顿）
     * new KalmanFilterUtil(2.0, 3.0);  // 中等参数，平衡跟踪精度
     * <p>
     * 参数调整原则
     * 运动规律性：越规整的运动，过程噪声设得越低
     * 信号质量：GPS信号越差，测量噪声设得越高
     * 实时性要求：对实时性要求高，适当提高过程噪声
     * 平滑度要求：对平滑度要求高，适当提高测量噪声
     *
     * @param processNoiseStd     过程噪声标准差
     * @param measurementNoiseStd 测量噪声标准差
     */
    public KalmanFilterUtil(double processNoiseStd, double measurementNoiseStd) {
        this.processNoiseStd = processNoiseStd;
        this.measurementNoiseStd = measurementNoiseStd;
        this.measurementNoise = new double[2][2];
        this.measurementNoise[0][0] = Math.pow(this.measurementNoiseStd, 2);
        this.measurementNoise[1][1] = Math.pow(this.measurementNoiseStd, 2);
        initialize();
    }

    public KalmanFilterUtil() {
        this(2.0, 8.0);
    }

    private void initialize() {
        stateMatrix = new double[][]{{0}, {0}, {0}, {0}};
        errorCovMatrix = new double[4][4];
        errorCovMatrix[0][0] = 10.0;  // x坐标误差（米）
        errorCovMatrix[1][1] = 10.0;  // y坐标误差（米）
        errorCovMatrix[2][2] = 5.0;   // vx速度误差（米/秒）
        errorCovMatrix[3][3] = 5.0;   // vy速度误差（米/秒）
    }

    public List<TrackPoint> filterTrack(List<TrackPoint> points) {
        List<TrackPoint> sortedPoints = points.stream()
                .filter(point -> point.getTime() != null)
                .sorted(Comparator.comparing(TrackPoint::getTime))
                .collect(Collectors.toList());

        if (sortedPoints.size() < 2) {
            return new ArrayList<>(sortedPoints);
        }

        List<TrackPoint> filteredPoints = new ArrayList<>();
        TrackPoint prevPoint = sortedPoints.get(0);

        // 初始化第一个点的状态
        Coordinate firstMercator = convertToMercator(prevPoint.getLon(), prevPoint.getLat());
        stateMatrix[0][0] = firstMercator.x;
        stateMatrix[1][0] = firstMercator.y;
        double firstSpeedMs = kmhToMs(prevPoint.getSpeed());
        double firstDirRadian = Math.toRadians(prevPoint.getDirection());

        // 速度分量计算（东向vx，北向vy）
        stateMatrix[2][0] = firstSpeedMs * Math.sin(firstDirRadian);
        stateMatrix[3][0] = firstSpeedMs * Math.cos(firstDirRadian);

        filteredPoints.add(prevPoint);

        for (int i = 1; i < sortedPoints.size(); i++) {
            TrackPoint currPoint = sortedPoints.get(i);
            double dt = Duration.between(prevPoint.getTime(), currPoint.getTime()).toMillis() / 1000.0;

            if (dt <= 0) {
                filteredPoints.add(currPoint);
                prevPoint = currPoint;
                continue;
            }

            predict(dt);

            Coordinate currMercator = convertToMercator(currPoint.getLon(), currPoint.getLat());
            update(currMercator.x, currMercator.y);

            // 生成滤波后的点
            Coordinate filteredWgs84 = convertToWgs84(stateMatrix[0][0], stateMatrix[1][0]);
            double filteredSpeedKmh = msToKmh(Math.sqrt(Math.pow(stateMatrix[2][0], 2) + Math.pow(stateMatrix[3][0], 2)));

            // 方向角反算（正北顺时针）
            double filteredDirRadian = Math.atan2(stateMatrix[2][0], stateMatrix[3][0]);
            double filteredDirDegree = Math.toDegrees(filteredDirRadian);
            filteredDirDegree = (filteredDirDegree + 360) % 360;

            filteredPoints.add(new TrackPoint(
                    filteredWgs84.x,
                    filteredWgs84.y,
                    currPoint.getTime(),
                    filteredSpeedKmh,
                    filteredDirDegree
            ));

            prevPoint = currPoint;
        }

        return filteredPoints;
    }

    private void predict(double dt) {
        double[][] transitionMatrix = new double[4][4];
        transitionMatrix[0][0] = 1;
        transitionMatrix[1][1] = 1;
        transitionMatrix[2][2] = 1;
        transitionMatrix[3][3] = 1;
        transitionMatrix[0][2] = dt;
        transitionMatrix[1][3] = dt;

        stateMatrix = multiplyMatrix(transitionMatrix, stateMatrix);

        double[][] F_T = transposeMatrix(transitionMatrix);
        double[][] P_temp = multiplyMatrix(transitionMatrix, errorCovMatrix);
        double[][] P_pred = multiplyMatrix(P_temp, F_T);
        errorCovMatrix = addMatrix(P_pred, getProcessNoiseMatrix(dt));
    }

    private void update(double measuredX, double measuredY) {
        double[][] measurementMatrix = {{1, 0, 0, 0}, {0, 1, 0, 0}};
        double[][] H_T = transposeMatrix(measurementMatrix);

        double[][] S_temp1 = multiplyMatrix(measurementMatrix, errorCovMatrix);
        double[][] S_temp2 = multiplyMatrix(S_temp1, H_T);
        double[][] S = addMatrix(S_temp2, measurementNoise);
        double[][] K_temp = multiplyMatrix(errorCovMatrix, H_T);
        double[][] K = multiplyMatrix(K_temp, invert2x2Matrix(S));

        double[][] z = {{measuredX}, {measuredY}};
        double[][] H_x = multiplyMatrix(measurementMatrix, stateMatrix);
        double[][] residual = subtractMatrix(z, H_x);
        stateMatrix = addMatrix(stateMatrix, multiplyMatrix(K, residual));

        double[][] I = {{1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}};
        double[][] K_H = multiplyMatrix(K, measurementMatrix);
        errorCovMatrix = multiplyMatrix(subtractMatrix(I, K_H), errorCovMatrix);
    }

    private double[][] getProcessNoiseMatrix(double dt) {
        double dt2 = dt * dt;
        double dt3 = dt2 * dt;
        double dt4 = dt3 * dt;
        double q = Math.pow(processNoiseStd, 2);

        double[][] processNoise = new double[4][4];
        processNoise[0][0] = dt4 / 4 * q;
        processNoise[0][2] = dt3 / 2 * q;
        processNoise[1][1] = dt4 / 4 * q;
        processNoise[1][3] = dt3 / 2 * q;
        processNoise[2][0] = dt3 / 2 * q;
        processNoise[2][2] = dt2 * q;
        processNoise[3][1] = dt3 / 2 * q;
        processNoise[3][3] = dt2 * q;

        return processNoise;
    }

    /**
     * WGS84经纬度 → 墨卡托坐标
     */
    private Coordinate convertToMercator(double lon, double lat) {
        try {
            Coordinate wgs84Coord = new Coordinate(lon, lat);
            Coordinate mercatorCoord = new Coordinate(); // 显式创建输出坐标对象
            JTS.transform(wgs84Coord, mercatorCoord, wgs84ToMercator); // 三参数正确调用方式
            return mercatorCoord;
        } catch (Exception e) {
            throw new RuntimeException("经纬度转墨卡托失败：" + lon + "," + lat, e);
        }
    }

    /**
     * 墨卡托坐标 → WGS84经纬度
     */
    private Coordinate convertToWgs84(double x, double y) {
        try {
            Coordinate mercatorCoord = new Coordinate(x, y);
            Coordinate wgs84Coord = new Coordinate(); // 显式创建输出坐标对象
            JTS.transform(mercatorCoord, wgs84Coord, mercatorToWgs84); // 三参数正确调用方式
            return wgs84Coord;
        } catch (Exception e) {
            throw new RuntimeException("墨卡托转经纬度失败：" + x + "," + y, e);
        }
    }

    private double kmhToMs(double kmh) {
        return kmh * 1000 / 3600;
    }

    private double msToKmh(double ms) {
        return ms * 3600 / 1000;
    }

    private double[][] multiplyMatrix(double[][] a, double[][] b) {
        int rowsA = a.length;
        int colsA = a[0].length;
        int colsB = b[0].length;

        double[][] result = new double[rowsA][colsB];
        for (int i = 0; i < rowsA; i++) {
            for (int k = 0; k < colsA; k++) {
                if (a[i][k] == 0) continue;
                for (int j = 0; j < colsB; j++) {
                    result[i][j] += a[i][k] * b[k][j];
                }
            }
        }
        return result;
    }

    private double[][] transposeMatrix(double[][] matrix) {
        int rows = matrix.length;
        int cols = matrix[0].length;

        double[][] transposed = new double[cols][rows];
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                transposed[j][i] = matrix[i][j];
            }
        }
        return transposed;
    }

    private double[][] addMatrix(double[][] a, double[][] b) {
        int rows = a.length;
        int cols = a[0].length;

        double[][] result = new double[rows][cols];
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                result[i][j] = a[i][j] + b[i][j];
            }
        }
        return result;
    }

    private double[][] subtractMatrix(double[][] a, double[][] b) {
        int rows = a.length;
        int cols = a[0].length;

        double[][] result = new double[rows][cols];
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                result[i][j] = a[i][j] - b[i][j];
            }
        }
        return result;
    }

    private double[][] invert2x2Matrix(double[][] matrix) {
        double det = matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];
        if (Math.abs(det) < 1e-10) {
            throw new ArithmeticException("矩阵不可逆，行列式接近0");
        }

        double invDet = 1.0 / det;
        double[][] inverse = new double[2][2];
        inverse[0][0] = matrix[1][1] * invDet;
        inverse[0][1] = -matrix[0][1] * invDet;
        inverse[1][0] = -matrix[1][0] * invDet;
        inverse[1][1] = matrix[0][0] * invDet;

        return inverse;
    }
}
