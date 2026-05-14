package sunyu.util.pojo;

import java.time.LocalDateTime;

/**
 * WGS84 坐标系 GPS 轨迹点实体类
 * <p>
 * 用于封装农机 GPS 设备采集的单个轨迹点的完整信息，包括地理坐标（经纬度）、
 * 定位时间、定位状态、作业状态和行驶速度等核心属性。
 * 该类是整个 GIS 工具链的基础数据单元，所有空间分析和几何计算均以此类为输入源。
 * </p>
 * <p>
 * <b>坐标系说明：</b>
 * 采用 WGS84（World Geodetic System 1984）地理坐标系（EPSG:4326），
 * 经度和纬度以度（°）为单位，是全球通用的 GPS 定位标准坐标系。
 * 该坐标系下的坐标可直接用于地图展示、地理编码和跨系统数据交换，
 * 但不适合直接进行距离计算和几何操作（需先转换为高斯投影平面坐标）。
 * </p>
 * <p>
 * <b>与 {@link GaussPoint} 的关系：</b>
 * <ul>
 *   <li>{@link Wgs84Point}：存储原始 GPS 经纬度坐标，用于数据输入、地图展示和序列化；</li>
 *   <li>{@link GaussPoint}：继承自本类，扩展了高斯投影平面坐标（X, Y），
 *       用于距离计算、几何缓冲、面积量算等空间分析操作。</li>
 * </ul>
 * </p>
 * <p>
 * <b>核心使用场景：</b>
 * <ul>
 *   <li>GPS 数据采集：接收并存储农机终端上报的原始定位数据；</li>
 *   <li>轨迹可视化：将经纬度坐标直接渲染到 Web 地图或 GIS 软件中；</li>
 *   <li>坐标转换输入：作为高斯投影转换的源数据，生成 {@link GaussPoint}；</li>
 *   <li>作业状态分析：根据 gpsStatus 和 jobStatus 判断点位有效性，过滤异常数据；</li>
 *   <li>速度统计：基于 speed 字段计算平均作业速度、识别超速或怠速。</li>
 * </ul>
 * </p>
 *
 * @author SunYu
 * @see GaussPoint
 * @see sunyu.util.GisUtil#toGaussPointList(java.util.List)
 * @since 1.0.0
 */
public class Wgs84Point {
    /**
     * GPS 定位时间
     * <p>
     * 记录该轨迹点的精确时间戳，由 GPS 卫星授时系统提供，精度通常为秒级或亚秒级。
     * 该字段是时序分析、轨迹排序、时间窗口拆分和作业时长统计的核心依据。
     * </p>
     * <p>
     * <b>时序约束：</b>同一批轨迹点应按 gpsTime 升序排列，
     * 相邻点的时间差反映作业连续性（正常作业间隔通常为 1~10 秒）。
     * </p>
     */
    private LocalDateTime gpsTime;
    /**
     * 经度（单位：度）
     * <p>
     * WGS84 坐标系下的经度值，以本初子午线为基准：
     * <ul>
     *   <li>东经为正值（0.0° 至 180.0°）；</li>
     *   <li>西经为负值（-180.0° 至 0.0°）。</li>
     * </ul>
     * </p>
     * <p>
     * <b>取值范围：</b>{@code -180.0 ≤ longitude ≤ 180.0}。
     * 超出此范围的值为非法坐标，应在数据校验阶段过滤。
     * </p>
     * <p>
     * <b>精度说明：</b>一般 GPS 设备精度为 3~10 米，对应经纬度小数点后 5~6 位；
     * 高精度 RTK 设备可达厘米级，对应小数点后 7~8 位。
     * </p>
     */
    private double longitude;
    /**
     * 纬度（单位：度）
     * <p>
     * WGS84 坐标系下的纬度值，以赤道为基准：
     * <ul>
     *   <li>北纬为正值（0.0° 至 90.0°）；</li>
     *   <li>南纬为负值（-90.0° 至 0.0°）。</li>
     * </ul>
     * </p>
     * <p>
     * <b>取值范围：</b>{@code -90.0 ≤ latitude ≤ 90.0}。
     * 超出此范围的值为非法坐标，应在数据校验阶段过滤。
     * </p>
     * <p>
     * <b>精度说明：</b>与经度一致，一般 GPS 设备精度为 3~10 米，对应小数点后 5~6 位。
     * </p>
     */
    private double latitude;
    /**
     * GPS 定位状态
     * <p>
     * 标识 GPS 模块当前是否成功获取卫星定位信号，用于判断坐标数据的可靠性。
     * </p>
     * <p>
     * <b>取值说明：</b>
     * <ul>
     *   <li>{@code 0} — 未知：定位状态未明确，通常为数据缺失或设备初始化阶段；</li>
     *   <li>{@code 1} — 已定位：GPS 成功获取卫星信号，坐标数据可信；</li>
     *   <li>{@code 2} — 未定位：GPS 未获取有效卫星信号（如室内、隧道、遮挡），坐标数据不可信，应过滤。</li>
     * </ul>
     * </p>
     * <p>
     * <b>处理建议：</b>在聚类和几何计算前，建议过滤掉 gpsStatus = 2 的点位，
     * 避免未定位点的漂移坐标污染地块轮廓。
     * </p>
     */
    private int gpsStatus;
    /**
     * 农机作业状态
     * <p>
     * 标识农机在该点位是否处于作业状态（如耕作、播种、收割等），
     * 由农机终端传感器（如机具升降传感器、PTO 转速传感器）自动判断并上报。
     * </p>
     * <p>
     * <b>取值说明：</b>
     * <ul>
     *   <li>{@code 0} — 未知：作业状态未明确，传感器数据缺失或设备未配置；</li>
     *   <li>{@code 1} — 作业中：农机机具处于工作状态，该点位属于有效作业轨迹；</li>
     *   <li>{@code 2} — 未作业：农机机具未工作（如空驶、转移、停车），该点位属于非作业轨迹。</li>
     * </ul>
     * </p>
     * <p>
     * 当 checkWorkingStatus 设置为 {@code true} 时，只有 jobStatus = 1 的点位会参与聚类和几何计算，
     * 有效排除道路行驶和田间转移对地块轮廓的干扰。
     * </p>
     */
    private int jobStatus;
    /**
     * GPS 速度（单位：千米/小时）
     * <p>
     * 记录农机在该点位时的瞬时行驶速度，由 GPS 模块根据多普勒频移计算得出。
     * 该字段可用于速度统计、作业效率评估和异常状态识别（如超速、怠速）。
     * </p>
     * <p>
     * <b>取值说明：</b>
     * <ul>
     *   <li>{@code null}：速度数据缺失，GPS 模块未上报速度信息；</li>
     *   <li>{@code ≥ 0}：有效速度值，单位 km/h。注意：静止时速度为 0，不代表数据异常。</li>
     * </ul>
     * </p>
     * <p>
     * <b>典型作业速度参考：</b>
     * <ul>
     *   <li>犁地：3~8 km/h</li>
     *   <li>播种：5~12 km/h</li>
     *   <li>收割：3~10 km/h</li>
     *   <li>道路转移：20~60 km/h</li>
     * </ul>
     * </p>
     */
    private Double speed;

    /**
     * 默认无参构造方法
     * <p>
     * 创建一个空的 Wgs84Point 实例，所有字段初始化为默认值（数值型为 0，引用型为 null）。
     * 通常用于反射创建、序列化反序列化或后续通过 setter 方法填充数据。
     * </p>
     */
    public Wgs84Point() {
    }

    /**
     * 经纬度构造方法
     * <p>
     * 创建一个仅包含地理坐标的 Wgs84Point 实例，适用于仅需位置信息而不需要时间、状态等附加属性的场景。
     * </p>
     *
     * @param longitude 经度（单位：度），范围 -180.0 至 180.0，东经为正，西经为负
     * @param latitude  纬度（单位：度），范围 -90.0 至 90.0，北纬为正，南纬为负
     */
    public Wgs84Point(double longitude, double latitude) {
        this.longitude = longitude;
        this.latitude = latitude;
    }

    /**
     * 时间+经纬度构造方法
     * <p>
     * 创建一个包含定位时间和地理坐标的 Wgs84Point 实例，
     * 适用于轨迹点的基础封装，支持时序分析和轨迹排序。
     * </p>
     *
     * @param gpsTime   GPS 定位时间，不可为 null
     * @param longitude 经度（单位：度），范围 -180.0 至 180.0，东经为正，西经为负
     * @param latitude  纬度（单位：度），范围 -90.0 至 90.0，北纬为正，南纬为负
     */
    public Wgs84Point(LocalDateTime gpsTime, double longitude, double latitude) {
        this.gpsTime = gpsTime;
        this.longitude = longitude;
        this.latitude = latitude;
    }

    /**
     * 时间+经纬度+定位状态构造方法
     * <p>
     * 创建一个包含定位时间、地理坐标和 GPS 定位状态的 Wgs84Point 实例，
     * 适用于需要判断坐标可靠性的场景（如过滤未定位点）。
     * </p>
     *
     * @param gpsTime   GPS 定位时间，不可为 null
     * @param longitude 经度（单位：度），范围 -180.0 至 180.0，东经为正，西经为负
     * @param latitude  纬度（单位：度），范围 -90.0 至 90.0，北纬为正，南纬为负
     * @param gpsStatus GPS 定位状态（0：未知；1：已定位；2：未定位）
     */
    public Wgs84Point(LocalDateTime gpsTime, double longitude, double latitude, int gpsStatus) {
        this.gpsTime = gpsTime;
        this.longitude = longitude;
        this.latitude = latitude;
        this.gpsStatus = gpsStatus;
    }

    /**
     * 全参构造方法（不含速度）
     * <p>
     * 创建一个包含定位时间、地理坐标、定位状态和作业状态的 Wgs84Point 实例，
     * 适用于完整的轨迹点封装，支持作业状态过滤和地块识别。
     * </p>
     *
     * @param gpsTime   GPS 定位时间，不可为 null
     * @param longitude 经度（单位：度），范围 -180.0 至 180.0，东经为正，西经为负
     * @param latitude  纬度（单位：度），范围 -90.0 至 90.0，北纬为正，南纬为负
     * @param gpsStatus GPS 定位状态（0：未知；1：已定位；2：未定位）
     * @param jobStatus 作业状态（0：未知；1：作业中；2：未作业）
     */
    public Wgs84Point(LocalDateTime gpsTime, double longitude, double latitude, int gpsStatus, int jobStatus) {
        this.gpsTime = gpsTime;
        this.longitude = longitude;
        this.latitude = latitude;
        this.gpsStatus = gpsStatus;
        this.jobStatus = jobStatus;
    }

    /**
     * 获取 GPS 定位时间
     *
     * @return GPS 定位时间（LocalDateTime）；若未设置则返回 null
     */
    public LocalDateTime getGpsTime() {
        return gpsTime;
    }

    /**
     * 设置 GPS 定位时间
     *
     * @param gpsTime GPS 定位时间，不可为 null
     */
    public void setGpsTime(LocalDateTime gpsTime) {
        this.gpsTime = gpsTime;
    }

    /**
     * 获取经度（单位：度）
     *
     * @return 经度值，范围 -180.0 至 180.0，东经为正，西经为负
     */
    public double getLongitude() {
        return longitude;
    }

    /**
     * 设置经度（单位：度）
     *
     * @param longitude 经度值，范围 -180.0 至 180.0，东经为正，西经为负
     */
    public void setLongitude(double longitude) {
        this.longitude = longitude;
    }

    /**
     * 获取纬度（单位：度）
     *
     * @return 纬度值，范围 -90.0 至 90.0，北纬为正，南纬为负
     */
    public double getLatitude() {
        return latitude;
    }

    /**
     * 设置纬度（单位：度）
     *
     * @param latitude 纬度值，范围 -90.0 至 90.0，北纬为正，南纬为负
     */
    public void setLatitude(double latitude) {
        this.latitude = latitude;
    }

    /**
     * 获取 GPS 定位状态
     *
     * @return GPS 定位状态（0：未知；1：已定位；2：未定位）
     */
    public int getGpsStatus() {
        return gpsStatus;
    }

    /**
     * 设置 GPS 定位状态
     *
     * @param gpsStatus GPS 定位状态（0：未知；1：已定位；2：未定位）
     */
    public void setGpsStatus(int gpsStatus) {
        this.gpsStatus = gpsStatus;
    }

    /**
     * 获取农机作业状态
     *
     * @return 作业状态（0：未知；1：作业中；2：未作业）
     */
    public int getJobStatus() {
        return jobStatus;
    }

    /**
     * 设置农机作业状态
     *
     * @param jobStatus 作业状态（0：未知；1：作业中；2：未作业）
     */
    public void setJobStatus(int jobStatus) {
        this.jobStatus = jobStatus;
    }

    /**
     * 获取 GPS 速度（单位：千米/小时）
     *
     * @return 速度值（km/h）；若未设置则返回 null
     */
    public Double getSpeed() {
        return speed;
    }

    /**
     * 设置 GPS 速度（单位：千米/小时）
     *
     * @param speed 速度值（km/h），必须大于等于 0；可为 null 表示速度数据缺失
     */
    public void setSpeed(Double speed) {
        this.speed = speed;
    }
}
