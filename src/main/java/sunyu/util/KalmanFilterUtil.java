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

public class KalmanFilterUtil {
    // 坐标系定义（不变）
    private static final CoordinateReferenceSystem WGS84_CRS = DefaultGeographicCRS.WGS84;
    private static final CoordinateReferenceSystem MERCATOR_CRS;
    private static MathTransform wgs84ToMercator;
    private static MathTransform mercatorToWgs84;

    // 卡尔曼滤波核心参数（核心修改：调整噪声参数）
    private double[][] stateMatrix;       // [x, y, vx, vy]
    private double[][] errorCovMatrix;    // 误差协方差矩阵 P
    private final double[][] measurementNoise;  // 测量噪声矩阵 R
    /**
     * 过程噪声标准差（增大到5.0）：让滤波器更信任预测（正常速度）
     * 越大→越信任测量，越小→越信任预测（当前设为5.0，平衡机动与平滑）
     */
    private double processNoiseStd = 5.0;
    /**
     * 测量噪声标准差（增大到50.0）：让滤波器少信任异常测量（偏移坐标）
     * 越大→越信任预测，越小→越信任测量（当前设为50.0，抗异常点）
     */
    private double measurementNoiseStd = 50.0;
    /**
     * 残差限制（米）：超过这个距离的测量偏差，视为异常，裁剪到阈值
     * 避免异常大的偏差过度影响速度计算（当前设为50米，根据需求调整）
     */
    private static final double RESIDUAL_LIMIT = 50.0;

    // 静态初始化坐标系（不变）
    static {
        try {
            MERCATOR_CRS = CRS.decode("EPSG:3857", true);
            wgs84ToMercator = CRS.findMathTransform(WGS84_CRS, MERCATOR_CRS, false);
            mercatorToWgs84 = CRS.findMathTransform(MERCATOR_CRS, WGS84_CRS, false);
            System.out.println("坐标系初始化成功：WGS84 → EPSG:3857");
        } catch (Exception e) {
            throw new RuntimeException("坐标系初始化失败", e);
        }
    }

    // 构造函数（同步更新默认参数）
    public KalmanFilterUtil(double processNoiseStd, double measurementNoiseStd) {
        this.processNoiseStd = processNoiseStd;
        this.measurementNoiseStd = measurementNoiseStd;
        this.measurementNoise = new double[2][2];
        this.measurementNoise[0][0] = Math.pow(this.measurementNoiseStd, 2);
        this.measurementNoise[1][1] = Math.pow(this.measurementNoiseStd, 2);
        initializeErrorCovMatrix();
    }

    public KalmanFilterUtil() {
        this(5.0, 50.0); // 默认使用调整后的参数
    }

    // 初始化误差协方差矩阵（不变）
    private void initializeErrorCovMatrix() {
        errorCovMatrix = new double[4][4];
        errorCovMatrix[0][0] = 10.0;  // x坐标误差（米）
        errorCovMatrix[1][1] = 10.0;  // y坐标误差（米）
        errorCovMatrix[2][2] = 8.0;   // vx速度误差（米/秒，匹配39.5km/h≈11m/s）
        errorCovMatrix[3][3] = 8.0;   // vy速度误差（米/秒）
    }

    // 对外暴露预处理方法（不变）
    public List<TrackPoint> getPreprocessedTrack(List<TrackPoint> originalPoints) {
        return preprocessTrackPoints(originalPoints);
    }

    // 预处理方法（不变，确保速度计算正确）
    private List<TrackPoint> preprocessTrackPoints(List<TrackPoint> originalPoints) {
        List<TrackPoint> sortedNewPoints = originalPoints.stream()
                .filter(point -> point.getTime() != null)
                .sorted(Comparator.comparing(TrackPoint::getTime))
                .map(point -> new TrackPoint(point.getLon(), point.getLat(), point.getTime()))
                .collect(Collectors.toList());

        if (sortedNewPoints.size() < 2) {
            return sortedNewPoints;
        }

        // 计算墨卡托坐标
        List<Coordinate> mercatorCoords = new ArrayList<>();
        for (TrackPoint point : sortedNewPoints) {
            Coordinate mercator = convertToMercator(point.getLon(), point.getLat());
            mercatorCoords.add(mercator);
            System.out.printf("预处理：第%d点 经纬度(%.6f,%.6f) → 墨卡托(%.1f,%.1f)%n",
                    sortedNewPoints.indexOf(point) + 1, point.getLon(), point.getLat(), mercator.x, mercator.y);
        }

        // 计算速度/方向角
        for (int i = 0; i < sortedNewPoints.size(); i++) {
            TrackPoint curr = sortedNewPoints.get(i);
            Coordinate currMercator = mercatorCoords.get(i);
            double speedKmh = 0.0;
            double direction = 0.0;

            if (i == 0) {
                TrackPoint next = sortedNewPoints.get(1);
                Coordinate nextMercator = mercatorCoords.get(1);
                double dt = Duration.between(curr.getTime(), next.getTime()).toMillis() / 1000.0;
                double distanceM = Math.sqrt(Math.pow(nextMercator.x - currMercator.x, 2) + Math.pow(nextMercator.y - currMercator.y, 2));
                speedKmh = (distanceM / dt) * 3.6;
                direction = calculateDirection(currMercator, nextMercator);
                System.out.printf("预处理：第1点 与第2点距离=%.1f米 → 速度=%.1f km/h%n", distanceM, speedKmh);
            } else {
                TrackPoint prev = sortedNewPoints.get(i - 1);
                Coordinate prevMercator = mercatorCoords.get(i - 1);
                double dt = Duration.between(prev.getTime(), curr.getTime()).toMillis() / 1000.0;
                double distanceM = Math.sqrt(Math.pow(currMercator.x - prevMercator.x, 2) + Math.pow(currMercator.y - prevMercator.y, 2));
                speedKmh = (distanceM / dt) * 3.6;
                direction = calculateDirection(prevMercator, currMercator);
                System.out.printf("预处理：第%d点 与第%d点距离=%.1f米 → 速度=%.1f km/h%n", i + 1, i, distanceM, speedKmh);
            }

            curr.setSpeed(speedKmh);
            curr.setDirection(direction);
        }

        return sortedNewPoints;
    }

    // 核心滤波方法（新增残差限制逻辑）
    public List<TrackPoint> filterTrack(List<TrackPoint> points) {
        if (points == null) {
            throw new IllegalArgumentException("输入轨迹点列表不能为null");
        }
        for (TrackPoint point : points) {
            if (point == null || point.getTime() == null) {
                throw new IllegalArgumentException("轨迹点及时间不能为null");
            }
        }

        // 使用预处理后的点
        List<TrackPoint> preprocessed = preprocessTrackPoints(points);
        if (preprocessed.size() < 2) {
            return new ArrayList<>(preprocessed);
        }

        List<TrackPoint> filtered = new ArrayList<>();
        TrackPoint first = preprocessed.get(0);
        Coordinate firstMercator = convertToMercator(first.getLon(), first.getLat());
        double firstSpeedMs = kmhToMs(first.getSpeed());
        double firstDirRadian = Math.toRadians(first.getDirection());

        // 初始化状态矩阵（用预处理的正常速度）
        stateMatrix = new double[4][1];
        stateMatrix[0][0] = firstMercator.x;
        stateMatrix[1][0] = firstMercator.y;
        stateMatrix[2][0] = firstSpeedMs * Math.sin(firstDirRadian); // vx≈11*sin(37.5°)≈6.7 m/s
        stateMatrix[3][0] = firstSpeedMs * Math.cos(firstDirRadian); // vy≈11*cos(37.5°)≈8.7 m/s

        // 添加第一个滤波点
        filtered.add(new TrackPoint(
                first.getLon(), first.getLat(), first.getTime(),
                first.getSpeed(), first.getDirection()
        ));

        // 处理后续点（核心优化：添加残差限制）
        for (int i = 1; i < preprocessed.size(); i++) {
            TrackPoint curr = preprocessed.get(i);
            TrackPoint prev = preprocessed.get(i - 1);
            double dt = Duration.between(prev.getTime(), curr.getTime()).toMillis() / 1000.0;
            if (dt <= 0.1) dt = 0.1;

            // 1. 预测步骤（不变）
            predict(dt);

            // 2. 更新步骤（新增：残差限制）
            Coordinate currMercator = convertToMercator(curr.getLon(), curr.getLat());
            updateWithResidualLimit(currMercator.x, currMercator.y); // 用新的更新方法

            // 3. 生成滤波点（不变）
            Coordinate filteredWgs84 = convertToWgs84(stateMatrix[0][0], stateMatrix[1][0]);
            double filteredSpeedMs = Math.sqrt(Math.pow(stateMatrix[2][0], 2) + Math.pow(stateMatrix[3][0], 2));
            // 可选：限制速度上限（避免极端值，如设为2倍正常速度）
            filteredSpeedMs = Math.min(filteredSpeedMs, firstSpeedMs * 2);
            double filteredDirDegree = (Math.toDegrees(Math.atan2(stateMatrix[2][0], stateMatrix[3][0])) + 360) % 360;

            filtered.add(new TrackPoint(
                    filteredWgs84.x, filteredWgs84.y, curr.getTime(),
                    msToKmh(filteredSpeedMs), filteredDirDegree
            ));
        }

        return filtered;
    }

    // ------------------------------ 新增：带残差限制的更新方法 ------------------------------

    /**
     * 带残差限制的更新步骤：超过阈值的位置偏差会被裁剪，避免异常影响速度
     */
    private void updateWithResidualLimit(double measuredX, double measuredY) {
        double[][] measurementMatrix = {{1, 0, 0, 0}, {0, 1, 0, 0}};
        double[][] H_T = transposeMatrix(measurementMatrix);

        // 计算卡尔曼增益K（不变）
        double[][] S_temp1 = multiplyMatrix(measurementMatrix, errorCovMatrix);
        double[][] S_temp2 = multiplyMatrix(S_temp1, H_T);
        double[][] S = addMatrix(S_temp2, measurementNoise);
        double[][] K_temp = multiplyMatrix(errorCovMatrix, H_T);
        double[][] K = multiplyMatrix(K_temp, invert2x2Matrix(S));

        // 计算残差（核心优化：限制残差大小）
        double[][] z = {{measuredX}, {measuredY}};
        double[][] H_x = multiplyMatrix(measurementMatrix, stateMatrix);
        double[][] residual = subtractMatrix(z, H_x);

        // 裁剪残差：x/y方向偏差超过RESIDUAL_LIMIT（50米）则限制在阈值内
        residual[0][0] = Math.max(Math.min(residual[0][0], RESIDUAL_LIMIT), -RESIDUAL_LIMIT);
        residual[1][0] = Math.max(Math.min(residual[1][0], RESIDUAL_LIMIT), -RESIDUAL_LIMIT);
        System.out.printf("更新：残差裁剪后 (x:%.1f米, y:%.1f米)%n", residual[0][0], residual[1][0]);

        // 用裁剪后的残差更新状态（避免速度异常）
        stateMatrix = addMatrix(stateMatrix, multiplyMatrix(K, residual));

        // 更新误差协方差（不变）
        double[][] I = {{1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}};
        double[][] K_H = multiplyMatrix(K, measurementMatrix);
        errorCovMatrix = multiplyMatrix(subtractMatrix(I, K_H), errorCovMatrix);
    }

    // ------------------------------ 其他方法（不变） ------------------------------
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

    private double calculateDirection(Coordinate prevMercator, Coordinate currMercator) {
        double deltaX = currMercator.x - prevMercator.x;
        double deltaY = currMercator.y - prevMercator.y;
        double radian = Math.atan2(deltaX, deltaY);
        return (Math.toDegrees(radian) + 360) % 360;
    }

    private Coordinate convertToMercator(double lon, double lat) {
        try {
            Coordinate wgs84 = new Coordinate(lon, lat);
            Coordinate mercator = new Coordinate();
            JTS.transform(wgs84, mercator, wgs84ToMercator);
            return mercator;
        } catch (Exception e) {
            throw new RuntimeException("经纬度(" + lon + "," + lat + ")转墨卡托失败", e);
        }
    }

    private Coordinate convertToWgs84(double x, double y) {
        try {
            Coordinate mercator = new Coordinate(x, y);
            Coordinate wgs84 = new Coordinate();
            JTS.transform(mercator, wgs84, mercatorToWgs84);
            return wgs84;
        } catch (Exception e) {
            throw new RuntimeException("墨卡托(" + x + "," + y + ")转经纬度失败", e);
        }
    }

    private double kmhToMs(double kmh) {
        return kmh / 3.6;
    }

    private double msToKmh(double ms) {
        return ms * 3.6;
    }

    // 矩阵运算方法（不变）
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
        if (Math.abs(det) < 1e-8) {
            throw new ArithmeticException("矩阵不可逆，行列式=" + det);
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