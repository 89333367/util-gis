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
    private final double processNoiseStd; // 过程噪声标准差
    private final double measurementNoiseStd; // 测量噪声标准差

    static {
        try {
            MERCATOR_CRS = CRS.decode("EPSG:3857");
            wgs84ToMercator = CRS.findMathTransform(WGS84_CRS, MERCATOR_CRS, true);
            mercatorToWgs84 = CRS.findMathTransform(MERCATOR_CRS, WGS84_CRS, true);
        } catch (Exception e) {
            throw new RuntimeException("坐标系初始化失败", e);
        }
    }

    public KalmanFilterUtil(double processNoiseStd, double measurementNoiseStd) {
        this.processNoiseStd = processNoiseStd;
        this.measurementNoiseStd = measurementNoiseStd;
        this.measurementNoise = new double[2][2];
        this.measurementNoise[0][0] = Math.pow(this.measurementNoiseStd, 2);
        this.measurementNoise[1][1] = Math.pow(this.measurementNoiseStd, 2);
        initialize();
    }

    public KalmanFilterUtil() {
        this(1.0, 5.0);
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
