package sunyu.util.test;

import cn.hutool.core.convert.Convert;
import cn.hutool.core.date.DateUtil;
import cn.hutool.core.date.LocalDateTimeUtil;
import cn.hutool.core.io.FileUtil;
import cn.hutool.core.io.resource.ResourceUtil;
import cn.hutool.core.util.ReUtil;
import cn.hutool.core.util.StrUtil;
import cn.hutool.db.Db;
import cn.hutool.db.DbUtil;
import cn.hutool.db.Entity;
import cn.hutool.log.Log;
import cn.hutool.log.LogFactory;
import cn.hutool.log.level.Level;
import com.zaxxer.hikari.HikariConfig;
import com.zaxxer.hikari.HikariDataSource;
import org.junit.jupiter.api.Test;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.MultiPolygon;
import sunyu.util.GisUtil;
import sunyu.util.pojo.*;

import javax.sql.DataSource;
import java.nio.charset.StandardCharsets;
import java.sql.SQLException;
import java.time.Duration;
import java.time.LocalDateTime;
import java.util.*;

public class TestUtilGis {
    Log log = LogFactory.get();
    GisUtil gisUtil = GisUtil.builder().build();
    String path = "D:/GitLab/util-gis/testFiles";
    private final String DID_REG = "[a-zA-Z0-9]+[ ]?[a-zA-Z0-9]+";

    private TreeMap<String, String> parseProtocolString(String protocolString) {
        if (StrUtil.isNotBlank(protocolString)) {
            String[] data = protocolString.split("\\$");
            // 1. 消息元数据间使用一个“$”为分隔符。
            // 2. 消息结构：消息前缀;序列号;终端 ID;命令标识;参数明细。固定 5 列。
            // 3. 前缀主要有 SUBMIT 和 REPORT 两种，SUBMIT 表示主动发送，REPORT 表示应答。
            // 4. 序列号主要用于下行指令的状态通知匹配，主动上传的消息序列号默认为 1。
            // 5. 国标终端的 VIN 用于终端 ID，其它终端 ID 为博创自定义 ID。
            // 6. 参数明细包含多个参数，单个参数以 KEY:VALUE 形式，多个参数以半角逗号分隔。
            // 7. 车辆登入示例： SUBMIT$1$LVBV4J0B2AJ063987$LOGIN$TIME:20150623120000,1001:1
            if (data.length == 5 && StrUtil.isNotBlank(data[2]) && ReUtil.isMatch(DID_REG, data[2])) {
                // 如果是内部协议的5段，并且，设备编号不为空
                TreeMap<String, String> params = new TreeMap<>();
                params.put("params0", data[0]);// 前缀,SUBMIT 和REPORT 两种，SUBMIT 表示主动发送，REPORT 表示应答
                params.put("params1", data[1]);// 消息id
                params.put("params2", data[2]);// 终端编号
                params.put("params3", data[3]);// 命令标识,PACKET,LINKSTATUS,TERMIN,TERMOUT,REALTIME,HISTORY,GETARG,SETARG,CONTROL,UPDATE,NOTIFY
                params.put("did", params.get("params2"));// 冗余key，便于检索
                String[] datas = data[4].split(",");// 数据部分
                for (String kv : datas) {// 拼装内部协议key数据
                    String[] keyValue = kv.split(":");
                    if (keyValue.length == 2) {// key和value都有，例如：key:value
                        String id = keyValue[0];// 内部协议数字key
                        String value = keyValue[1];// 网关转换后的值
                        if (StrUtil.isNotBlank(id) && StrUtil.isNotBlank(value)) {
                            params.put(id, value);// 使用内部协议的数字key
                        }
                    }
                }
                return params;
            }
        }
        return null;
    }

    private DataSource getMySqlDatasource() {
        HikariConfig config = new HikariConfig();
        config.setJdbcUrl("jdbc:mysql://172.16.1.59:3306/farm?useUnicode=true&characterEncoding=UTF-8&serverTimezone=UTC&zeroDateTimeBehavior=convertToNull&useInformationSchema=true&useSSL=false&allowMultiQueries=true");
        config.setUsername("dev");
        config.setPassword("uml-tech");
        config.setMinimumIdle(0);
        config.setMaximumPoolSize(10);
        DataSource ds = new HikariDataSource(config);
        return ds;
    }

    private Db getMysqlDb() {
        DataSource ds = getMySqlDatasource();
        DbUtil.setShowSqlGlobal(true, false, true, Level.DEBUG);
        Db db = Db.use(ds);
        return db;
    }

    private HikariDataSource getTdengineDatasource() {
        HikariConfig config = new HikariConfig();
        config.setDriverClassName("com.taosdata.jdbc.ws.WebSocketDriver");
        config.setJdbcUrl("jdbc:TAOS-WS://172.16.1.173:16041/?httpConnectTimeout=60000&messageWaitTimeout=60000");
        config.setUsername("root");
        config.setPassword("taosdata");
        config.setMinimumIdle(0);
        config.setMaximumPoolSize(10);
        HikariDataSource ds = new HikariDataSource(config);
        return ds;
    }

    private Db getTdengineDb() {
        DataSource ds = getTdengineDatasource();
        DbUtil.setShowSqlGlobal(true, false, true, Level.DEBUG);
        Db db = Db.use(ds);
        return db;
    }

    void 生成数据文件(String did, String startTime, String endTime) {
        if (!FileUtil.exist(path + StrUtil.format("/{}_{}_{}_trace.txt", did, startTime, endTime))) {
            Db tdengineDb = getTdengineDb();
            String jobStartTime = DateUtil.parse(startTime, "yyyyMMddHHmmss").toString("yyyy-MM-dd HH:mm:ss");
            String jobEndTime = DateUtil.parse(endTime, "yyyyMMddHHmmss").toString("yyyy-MM-dd HH:mm:ss");
            String tdSql = StrUtil.format("select protocol from frequent.d_p where did='{}' and _rowts>='{}' and _rowts<='{}'", did, jobStartTime, jobEndTime);
            log.debug("{}", tdSql);
            List<Entity> rows = null;
            try {
                rows = tdengineDb.query(tdSql);
            } catch (SQLException e) {
                throw new RuntimeException(e);
            }
            List<String> traceList = new ArrayList<>();
            List<String> protocolList = new ArrayList<>();
            for (Entity row : rows) {
                String protocolStr = new String(row.getBytes("protocol"), StandardCharsets.UTF_8);
                Map<String, String> protocol = parseProtocolString(protocolStr);
                if (!protocol.containsKey("3014") || !protocol.containsKey("2602") || !protocol.containsKey("2603")) {
                    continue;
                }
                protocolList.add(protocolStr);//只要有定位时间和经纬度信息就添加到协议列表
                if (protocol.containsKey("2601") && !protocol.get("2601").equals("0")) {// 定位状态,0已定位，1未定位
                    continue;
                }
                // 经纬度不能为0（无效坐标）
                if (Convert.toDouble(protocol.get("2602")) == 0.0 && Convert.toDouble(protocol.get("2603")) == 0.0) {
                    continue;
                }
                // 经纬度必须在合理范围内
                if (Convert.toDouble(protocol.get("2602")) < -180.0 || Convert.toDouble(protocol.get("2602")) > 180.0 || Convert.toDouble(protocol.get("2603")) < -90.0 || Convert.toDouble(protocol.get("2603")) > 90.0) {
                    continue;
                }
                // {定位时间yyyyMMddHHmmss},{经度},{纬度}
                traceList.add(StrUtil.format("{},{},{}", protocol.get("3014"), protocol.get("2602"), protocol.get("2603")));
            }
            FileUtil.writeUtf8Lines(traceList, path + StrUtil.format("/{}_{}_{}_trace.txt", did, startTime, endTime));
            FileUtil.writeUtf8Lines(protocolList, path + StrUtil.format("/{}_{}_{}_protocol.txt", did, startTime, endTime));
        }
    }

    void 测试拆分数据(String did, String startTime, String endTime, double jobWidth) {
        String fileName = path + StrUtil.format("/{}_{}_{}_protocol.txt", did, startTime, endTime);
        if (!FileUtil.exist(fileName)) {
            return;
        }
        List<Wgs84Point> l = new ArrayList<>();
        for (String protocolStr : FileUtil.readUtf8Lines(fileName)) {
            Map<String, String> protocol = parseProtocolString(protocolStr);
            Wgs84Point wgs84Point = new Wgs84Point();
            wgs84Point.setGpsTime(LocalDateTimeUtil.parse(protocol.get("3014"), "yyyyMMddHHmmss"));
            wgs84Point.setLongitude(Double.parseDouble(protocol.get("2602")));
            wgs84Point.setLatitude(Double.parseDouble(protocol.get("2603")));
            if (protocol.containsKey("2601")) {// 定位状态,0已定位，1未定位
                if (protocol.get("2601").equals("0")) {
                    wgs84Point.setGpsStatus(1);
                } else {
                    wgs84Point.setGpsStatus(2);
                }
            }
            if (protocol.containsKey("3020")) {// 终端ACC状态,0关闭，1开启
                if (Convert.toStr(protocol.get("3020")).equals("0")) {
                    wgs84Point.setJobStatus(2);//ACC关闭，认为是没有作业
                }
            }
            if (protocol.containsKey("4031")) {// 作业标识,1作业,0非作业,2暂停
                if (!Convert.toStr(protocol.get("4031")).equals("1")) {
                    wgs84Point.setJobStatus(2);//作业标识不是1，认为是没有作业
                }
            }
            l.add(wgs84Point);
        }

        String partsFile = StrUtil.format(path + "/{}_{}_{}_parts.txt", did, startTime, endTime);
        List<String> partsInfo = new ArrayList<>();
        SplitResult splitResult = gisUtil.splitRoad(l, jobWidth);
        partsInfo.add(StrUtil.format("作业总幅宽（米）: {}", splitResult.getWorkingWidth()));
        partsInfo.add(StrUtil.format("总WKT: {}", splitResult.getWkt()));
        partsInfo.add(StrUtil.format("作业总面积（亩）: {}", splitResult.getMu()));
        partsInfo.add(StrUtil.format("作业时间范围: {} - {}", splitResult.getStartTime(), splitResult.getEndTime()));
        partsInfo.add(StrUtil.format("最小有效时间间隔（秒）: {}", splitResult.getMinEffectiveInterval()));
        partsInfo.add(StrUtil.format("合并后有 {} 个地块", splitResult.getGaussGeometry().getNumGeometries()));
        partsInfo.add(StrUtil.format("拆分后有 {} 个地块", splitResult.getParts().size()));
        partsInfo.add("\n");
        int partIndex = 1;
        for (SplitPart splitPart : splitResult.getParts()) {
            List<String> partInfo = new ArrayList<>();
            partInfo.add(StrUtil.format("地块 {}:", partIndex++));
            partInfo.add(StrUtil.format("子WKT: {}", splitPart.getWkt()));
            partInfo.add(StrUtil.format("作业面积（亩）: {}", splitPart.getMu()));
            partInfo.add(StrUtil.format("作业时间范围: {} - {}", splitPart.getStartTime(), splitPart.getEndTime()));
            partInfo.add(StrUtil.format("最小有效时间间隔（秒）: {}", splitPart.getMinEffectiveInterval()));
            if (splitPart.getGaussGeometry() instanceof MultiPolygon) {
                partInfo.add(StrUtil.format("有 {} 个子地块", splitPart.getGaussGeometry().getNumGeometries()));
            }
            partsInfo.add(StrUtil.join("\n", partInfo) + "\n");
        }
        FileUtil.writeUtf8Lines(partsInfo, partsFile);

        fileName = path + StrUtil.format("/{}_{}_{}_gauss.txt", did, startTime, endTime);
        l = gisUtil.filterWgs84Points(l);
        List<GaussPoint> gaussPointList = gisUtil.toGaussPointList(l);
        List<String> gaussXyList = new ArrayList<>();
        for (GaussPoint gaussPoint : gaussPointList) {
            gaussXyList.add(StrUtil.format("{},{}", gaussPoint.getGaussX(), gaussPoint.getGaussY()));
        }
        FileUtil.writeUtf8Lines(gaussXyList, fileName);
    }

    void 生成HTML(String did, String startTime, String endTime) {
        String html = ResourceUtil.readUtf8Str("showGeometrysTemplate.html");

        String fileName = path + StrUtil.format("/{}_{}_{}_trace.txt", did, startTime, endTime);
        if (!FileUtil.exist(fileName)) {
            return;
        }
        String trace = FileUtil.readUtf8String(fileName);
        html = StrUtil.replace(html, "${trace}", trace);

        List<String> outlineList = new ArrayList<>();
        for (String line : FileUtil.readUtf8Lines(path + StrUtil.format("/{}_{}_{}_parts.txt", did, startTime, endTime))) {
            if (line.startsWith("子WKT")) {
                outlineList.add(line.replace("子WKT: ", ""));
            }
        }
        if (!outlineList.isEmpty()) {
            html = StrUtil.replace(html, "${outline}", StrUtil.join("\n", outlineList));
        } else {
            String outline = FileUtil.readUtf8Lines(path + StrUtil.format("/{}_{}_{}_parts.txt", did, startTime, endTime)).get(1);
            html = StrUtil.replace(html, "${outline}", outline.replace("总WKT: ", ""));
        }

        FileUtil.writeUtf8String(html, path + StrUtil.format("/{}_{}_{}.html", did, startTime, endTime));
    }

    void 查看速度(String did, String startTime, String endTime, LocalDateTime aTime, LocalDateTime bTime) {
        String fileName = path + StrUtil.format("/{}_{}_{}_trace.txt", did, startTime, endTime);
        if (!FileUtil.exist(fileName)) {
            return;
        }
        List<String> lines = FileUtil.readUtf8Lines(fileName);
        List<Wgs84Point> points = new ArrayList<>();
        for (String line : lines) {
            String[] info = line.split(",");
            LocalDateTime gpsTime = LocalDateTimeUtil.parse(info[0], "yyyyMMddHHmmss");
            double longitude = Double.parseDouble(info[1]);
            double latitude = Double.parseDouble(info[2]);
            if (aTime.isBefore(gpsTime) && gpsTime.isBefore(bTime)) {
                points.add(new Wgs84Point(gpsTime, longitude, latitude));
            }
        }
        // 打印 最小速度、最大速度、平均速度、速度中位数
        double[] speeds = new double[points.size() - 1];
        for (int i = 0; i < speeds.length; i++) {
            Wgs84Point p1 = points.get(i);
            Wgs84Point p2 = points.get(i + 1);
            double dist = gisUtil.haversine(p1, p2);          // 米
            double dt = Duration.between(p1.getGpsTime(), p2.getGpsTime()).getSeconds(); // 秒
            speeds[i] = dt == 0 ? 0 : (dist / 1000.0) / (dt / 3600.0); // km/h
        }
        Arrays.sort(speeds);
        double min = speeds.length > 0 ? speeds[0] : 0;
        double max = speeds.length > 0 ? speeds[speeds.length - 1] : 0;
        double avg = Arrays.stream(speeds).average().orElse(0);
        double median = speeds.length % 2 == 0
                ? (speeds[speeds.length / 2 - 1] + speeds[speeds.length / 2]) / 2
                : speeds[speeds.length / 2];
        log.info("速度统计  min={} km/h, max={} km/h, avg={} km/h, median={} km/h",
                String.format("%.2f", min),
                String.format("%.2f", max),
                String.format("%.2f", avg),
                String.format("%.2f", median));
    }

    @Test
    void 测试内存与cpu占用() {
        String did = "EC71BT2406060220";
        String startTime = "20251102154200";
        String endTime = "20251102172202";
        double jobWidth = 1.75;
        String fileName = path + StrUtil.format("/{}_{}_{}_protocol.txt", did, startTime, endTime);
        if (!FileUtil.exist(fileName)) {
            return;
        }
        List<Wgs84Point> l = new ArrayList<>();
        for (String protocolStr : FileUtil.readUtf8Lines(fileName)) {
            Map<String, String> protocol = parseProtocolString(protocolStr);
            Wgs84Point wgs84Point = new Wgs84Point();
            wgs84Point.setGpsTime(LocalDateTimeUtil.parse(protocol.get("3014"), "yyyyMMddHHmmss"));
            wgs84Point.setLongitude(Double.parseDouble(protocol.get("2602")));
            wgs84Point.setLatitude(Double.parseDouble(protocol.get("2603")));
            if (protocol.containsKey("2601")) {// 定位状态,0已定位，1未定位
                if (protocol.get("2601").equals("0")) {
                    wgs84Point.setGpsStatus(1);
                } else {
                    wgs84Point.setGpsStatus(2);
                }
            }
            if (protocol.containsKey("3020")) {// 终端ACC状态,0关闭，1开启
                if (Convert.toStr(protocol.get("3020")).equals("0")) {
                    wgs84Point.setJobStatus(2);//ACC关闭，认为是没有作业
                }
            }
            if (protocol.containsKey("4031")) {// 作业标识,1作业,0非作业,2暂停
                if (!Convert.toStr(protocol.get("4031")).equals("1")) {
                    wgs84Point.setJobStatus(2);//作业标识不是1，认为是没有作业
                }
            }
            l.add(wgs84Point);
        }
        for (int i = 0; i < 100; i++) {
            gisUtil.splitRoad(l, jobWidth);
        }
    }

    @Test
    void 测试批量投影转换() {
        String fileName = path + "/EC73BD2509061335_20251104100606_20251104101419_protocol.txt";
        List<Wgs84Point> l = new ArrayList<>();
        for (String protocolStr : FileUtil.readUtf8Lines(fileName)) {
            Map<String, String> protocol = parseProtocolString(protocolStr);
            Wgs84Point wgs84Point = new Wgs84Point();
            wgs84Point.setGpsTime(LocalDateTimeUtil.parse(protocol.get("3014"), "yyyyMMddHHmmss"));
            wgs84Point.setLongitude(Double.parseDouble(protocol.get("2602")));
            wgs84Point.setLatitude(Double.parseDouble(protocol.get("2603")));
            if (protocol.containsKey("2601")) {// 定位状态,0已定位，1未定位
                if (protocol.get("2601").equals("0")) {
                    wgs84Point.setGpsStatus(1);
                } else {
                    wgs84Point.setGpsStatus(2);
                }
            }
            if (protocol.containsKey("3020")) {// 终端ACC状态,0关闭，1开启
                if (Convert.toStr(protocol.get("3020")).equals("0")) {
                    wgs84Point.setJobStatus(2);//ACC关闭，认为是没有作业
                }
            }
            if (protocol.containsKey("4031")) {// 作业标识,1作业,0非作业,2暂停
                if (!Convert.toStr(protocol.get("4031")).equals("1")) {
                    wgs84Point.setJobStatus(2);//作业标识不是1，认为是没有作业
                }
            }
            l.add(wgs84Point);
        }
        List<Wgs84Point> wgs84Points = gisUtil.filterWgs84Points(l);
        log.debug("过滤后的wgs84点");
        for (Wgs84Point wgs84Point : wgs84Points) {
            log.debug("{} {} {}", wgs84Point.getGpsTime(), wgs84Point.getLongitude(), wgs84Point.getLatitude());
        }
        List<GaussPoint> gaussPoints = gisUtil.toGaussPointList(wgs84Points);
        log.debug("转换后的高斯投影点");
        for (GaussPoint gaussPoint : gaussPoints) {
            log.debug("{} {} {}", gaussPoint.getGpsTime(), gaussPoint.getGaussX(), gaussPoint.getGaussY());
        }
        List<Wgs84Point> wgs84Points1 = gisUtil.toWgs84PointList(gaussPoints);
        log.debug("转换后的wgs84点");
        for (Wgs84Point wgs84Point : wgs84Points1) {
            log.debug("{} {} {}", wgs84Point.getGpsTime(), wgs84Point.getLongitude(), wgs84Point.getLatitude());
        }
        List<Wgs84Point> closestPoints = gisUtil.findClosestPointList(wgs84Points1, wgs84Points);
        log.debug("容差点");
        for (Wgs84Point closestPoint : closestPoints) {
            log.debug("{} {} {}", closestPoint.getGpsTime(), closestPoint.getLongitude(), closestPoint.getLatitude());
        }
    }

    @Test
    void 测试单个点的投影转换1() {
        Wgs84Point p = new Wgs84Point(107.79342153, 22.54082875);
        log.debug("原始wgs84点");
        log.info("{} {} {}", p.getGpsTime(), p.getLongitude(), p.getLatitude());
        GaussPoint gaussPoint = gisUtil.toGaussPointList(new ArrayList<Wgs84Point>() {{
            add(p);
        }}).get(0);
        log.debug("转换后的高斯投影点");
        log.info("{} {} {}", gaussPoint.getGpsTime(), gaussPoint.getGaussX(), gaussPoint.getGaussY());
        Wgs84Point p1 = gisUtil.toWgs84PointList(new ArrayList<GaussPoint>() {{
            add(gaussPoint);
        }}).get(0);
        log.debug("转换后的wgs84点");
        log.info("{} {} {}", p1.getGpsTime(), p1.getLongitude(), p1.getLatitude());
        Wgs84Point closestPoint = gisUtil.findClosestPoint(p1, new ArrayList<Wgs84Point>() {{
            add(p);
        }}, 0.1);
        log.debug("找到最接近的wgs84点");
        log.info("{} {} {}", closestPoint.getGpsTime(), closestPoint.getLongitude(), closestPoint.getLatitude());
    }

    @Test
    void 测试点是否在圆中() {
        Wgs84Point p = new Wgs84Point(100.401807, 23.443696);
        Wgs84Point center = new Wgs84Point(100.27786, 23.60424);
        double radius = 1000.0;//米
        boolean isIn = gisUtil.inCircle(p, center, radius);
        log.info("点是否在圆中: {}", isIn);
    }

    @Test
    void 测试点是否在矩形中() {
        Wgs84Point p = new Wgs84Point(116.55470301, 40.21296700);
        Wgs84Point p1 = new Wgs84Point(116.55560000, 40.21296700); // 向东偏移约100米
        Wgs84Point p2 = new Wgs84Point(116.55560000, 40.21364248); // 向北偏移约100米
        boolean isIn = gisUtil.inRectangle(p, p1, p2);
        log.info("点是否在矩形中: {}", isIn);
    }

    @Test
    void 测试点是否在多边形中() {
        Geometry geom = gisUtil.toWgs84Geometry("POLYGON((116.55470301 40.21296700, 116.55560000 40.21296700, 116.55560000 40.21364248, 116.55470301 40.21364248, 116.55470301 40.21296700))");

        // 测试多边形内部的点
        Wgs84Point p1 = new Wgs84Point(116.55515000, 40.21330000);
        boolean isIn1 = gisUtil.inGeometry(p1, geom);
        log.info("内部点[116.55515000, 40.21330000]是否在多边形中: {}", isIn1);

        // 测试多边形顶点（边界点）
        Wgs84Point p2 = new Wgs84Point(116.55560000, 40.21364248);
        boolean isIn2 = gisUtil.inGeometry(p2, geom);
        log.info("顶点[116.55560000, 40.21364248]是否在多边形中: {}", isIn2);

        // 测试多边形外部的点
        Wgs84Point p3 = new Wgs84Point(116.55600000, 40.21400000);
        boolean isIn3 = gisUtil.inGeometry(p3, geom);
        log.info("外部点[116.55600000, 40.21400000]是否在多边形中: {}", isIn3);
    }

    @Test
    void 测试两点距离() {
        Wgs84Point p1 = new Wgs84Point(100.401807, 23.443696);
        Wgs84Point p2 = new Wgs84Point(100.27786, 23.60424);
        double distance = gisUtil.haversine(p1, p2);
        log.info("{} 米", distance);

        p1 = new Wgs84Point(116.55470301, 40.21296700);
        p2 = new Wgs84Point(116.55473883, 40.21364248);
        distance = gisUtil.haversine(p1, p2);
        log.info("{} 米", distance);
    }

    @Test
    void 测试1秒间隔001() {
        // 一块地，有镂空
        String did = "EC71BT2406060220";
        String startTime = "20251102154200";
        String endTime = "20251102172202";
        double jobWidth = 1.75;
        生成数据文件(did, startTime, endTime);
        测试拆分数据(did, startTime, endTime, jobWidth);
        生成HTML(did, startTime, endTime);
    }

    @Test
    void 测试1秒间隔002() {
        // 三块地，但有两块地之间的路的距离过小，所以合并成两块地
        String did = "EC71BT2406060220";
        String startTime = "20251102130028";
        String endTime = "20251102153804";
        double jobWidth = 1.75;
        生成数据文件(did, startTime, endTime);
        测试拆分数据(did, startTime, endTime, jobWidth);
        生成HTML(did, startTime, endTime);
    }

    @Test
    void 测试1秒间隔003() {
        // 一块地，轨迹很乱
        String did = "EC73BD2509060398";
        String startTime = "20251103130852";
        String endTime = "20251103151309";
        double jobWidth = 1.0;
        生成数据文件(did, startTime, endTime);
        测试拆分数据(did, startTime, endTime, jobWidth);
        生成HTML(did, startTime, endTime);
    }

    @Test
    void 测试1秒间隔004() {
        // 识别出了三块地，有镂空，有大石头
        String did = "EC71BT2406060220";
        String startTime = "20251103102528";
        String endTime = "20251103150242";
        double jobWidth = 1.75;
        生成数据文件(did, startTime, endTime);
        测试拆分数据(did, startTime, endTime, jobWidth);
        生成HTML(did, startTime, endTime);
    }

    @Test
    void 测试1秒间隔005() {
        // 广西测试地块，有作业开关状态上报
        String did = "EC73BD2506050018";
        String startTime = "20251104090717";
        String endTime = "20251104092257";
        double jobWidth = 2.8;
        生成数据文件(did, startTime, endTime);
        测试拆分数据(did, startTime, endTime, jobWidth);
        生成HTML(did, startTime, endTime);
    }

    @Test
    void 测试1秒间隔006() {
        // 广西测试地块，有作业开关状态上报
        String did = "EC73BD2509061335";
        String startTime = "20251104100606";
        String endTime = "20251104101419";
        double jobWidth = 2.8;
        生成数据文件(did, startTime, endTime);
        测试拆分数据(did, startTime, endTime, jobWidth);
        生成HTML(did, startTime, endTime);
    }

    @Test
    void 计算重复亩数0018_1335() {
        String wkt1 = FileUtil.readUtf8Lines(path + "/EC73BD2506050018_20251104090717_20251104092257_parts.txt").get(1).replace("总WKT: ", "");
        String wkt2 = FileUtil.readUtf8Lines(path + "/EC73BD2509061335_20251104100606_20251104101419_parts.txt").get(1).replace("总WKT: ", "");
        WktIntersectionResult r = gisUtil.intersection(wkt1, wkt2);
        log.info("相交轮廓WKT: {}", r.getWkt());
        log.info("相交面积：{} 亩", r.getMu());
    }

    @Test
    void 测试1秒间隔007() {
        // 两块地，有很长的路，有一块地蝴蝶形状
        String did = "EC73BD2509060248";
        String yyyyMMdd = "20251023";
        String startTime = yyyyMMdd + "000000";
        String endTime = yyyyMMdd + "235959";
        double jobWidth = 2.5;
        生成数据文件(did, startTime, endTime);
        测试拆分数据(did, startTime, endTime, jobWidth);
        生成HTML(did, startTime, endTime);
    }

    @Test
    void 测试1秒间隔008() {
        // 一天干了好多块地
        String did = "EC73BD2508220055";
        String yyyyMMdd = "20251013";
        String startTime = yyyyMMdd + "000000";
        String endTime = yyyyMMdd + "235959";
        double jobWidth = 2.6;
        生成数据文件(did, startTime, endTime);
        测试拆分数据(did, startTime, endTime, jobWidth);
        生成HTML(did, startTime, endTime);
    }

    @Test
    void 测试1秒间隔009() {
        // 有两块地，非常多的路
        String did = "NJ4GBQSAX0000687";
        String yyyyMMdd = "20250503";
        String startTime = yyyyMMdd + "000000";
        String endTime = yyyyMMdd + "235959";
        double jobWidth = 2.7;
        生成数据文件(did, startTime, endTime);
        测试拆分数据(did, startTime, endTime, jobWidth);
        生成HTML(did, startTime, endTime);
    }

    @Test
    void 测试1秒间隔010() {
        // 不会形成作业地块，不是正常的作业轨迹
        String did = "EC73BD2509061335";
        String startTime = "20251029101603";
        String endTime = "20251029102521";
        double jobWidth = 2;
        生成数据文件(did, startTime, endTime);
        测试拆分数据(did, startTime, endTime, jobWidth);
        生成HTML(did, startTime, endTime);
    }

    @Test
    void 测试1秒间隔011() {
        // 高密度地块，三块地，路稍微有一点点粘连
        String did = "NJ4GNBSAX0000693";
        String startTime = "20250505080044";
        String endTime = "20250505173658";
        double jobWidth = 2.5;
        生成数据文件(did, startTime, endTime);
        测试拆分数据(did, startTime, endTime, jobWidth);
        生成HTML(did, startTime, endTime);
    }

    @Test
    void 测试1秒间隔012() {
        // 高密度地块，三块地，路有一些粘连，地中有一点断开，地中有一点缝隙
        String did = "NJ4GNBSAX0000693";
        String startTime = "20250507073041";
        String endTime = "20250507162457";
        double jobWidth = 2.5;
        生成数据文件(did, startTime, endTime);
        测试拆分数据(did, startTime, endTime, jobWidth);
        生成HTML(did, startTime, endTime);
    }

    @Test
    void 测试1秒间隔013() {
        // 一块地，有很长的路
        String did = "NJ4GNBSAX0000693";
        String startTime = "20250509092721";
        String endTime = "20250509111620";
        double jobWidth = 2.5;
        生成数据文件(did, startTime, endTime);
        测试拆分数据(did, startTime, endTime, jobWidth);
        生成HTML(did, startTime, endTime);
    }

    @Test
    void 测试1秒间隔014() {
        // 时间有交叉
        String did = "EC71BT2404140062";
        String yyyyMMdd = "20240426";
        String startTime = yyyyMMdd + "000000";
        String endTime = yyyyMMdd + "235959";
        double jobWidth = 3.5;
        生成数据文件(did, startTime, endTime);
        测试拆分数据(did, startTime, endTime, jobWidth);
        生成HTML(did, startTime, endTime);
    }

    @Test
    void 测试1秒间隔015() {
        // 骗补
        String did = "NJ4GNBZAX0000160";
        String startTime = "20250603144925";
        String endTime = "20250603155017";
        double jobWidth = 2.4;
        生成数据文件(did, startTime, endTime);
        测试拆分数据(did, startTime, endTime, jobWidth);
        生成HTML(did, startTime, endTime);
    }

    @Test
    void 测试1秒间隔016() {
        String did = "NJ4GNBZAX0000273";
        String yyyyMMdd = "20250531";
        String startTime = yyyyMMdd + "000000";
        String endTime = yyyyMMdd + "235959";
        double jobWidth = 2.3;
        生成数据文件(did, startTime, endTime);
        测试拆分数据(did, startTime, endTime, jobWidth);
        生成HTML(did, startTime, endTime);
    }

    @Test
    void 测试1秒间隔017() {
        String did = "EC71BT2404140062";
        String yyyyMMdd = "20240427";
        String startTime = yyyyMMdd + "000000";
        String endTime = yyyyMMdd + "235959";
        double jobWidth = 2.1;
        生成数据文件(did, startTime, endTime);
        测试拆分数据(did, startTime, endTime, jobWidth);
        生成HTML(did, startTime, endTime);
    }

    @Test
    void 测试1秒间隔018() {
        String did = "EC71BT2404140062";
        String yyyyMMdd = "20240427";
        String startTime = yyyyMMdd + "000000";
        String endTime = yyyyMMdd + "235959";
        double jobWidth = 2.1;
        生成数据文件(did, startTime, endTime);
        测试拆分数据(did, startTime, endTime, jobWidth);
        生成HTML(did, startTime, endTime);
    }

    @Test
    void 测试1秒间隔019() {
        String did = "EC71BT2404140062";
        String yyyyMMdd = "20240428";
        String startTime = yyyyMMdd + "000000";
        String endTime = yyyyMMdd + "235959";
        double jobWidth = 3.9;
        生成数据文件(did, startTime, endTime);
        测试拆分数据(did, startTime, endTime, jobWidth);
        生成HTML(did, startTime, endTime);
    }

    @Test
    void 测试1秒间隔020() {
        String did = "EC71BT2404140062";
        String yyyyMMdd = "20240429";
        String startTime = yyyyMMdd + "000000";
        String endTime = yyyyMMdd + "235959";
        double jobWidth = 3.9;
        生成数据文件(did, startTime, endTime);
        测试拆分数据(did, startTime, endTime, jobWidth);
        生成HTML(did, startTime, endTime);
    }

    @Test
    void 测试1秒间隔021() {
        String did = "JFT3352503S00207";
        String yyyyMMdd = "20251024";
        String startTime = yyyyMMdd + "000000";
        String endTime = yyyyMMdd + "235959";
        double jobWidth = 3.5;
        生成数据文件(did, startTime, endTime);
        测试拆分数据(did, startTime, endTime, jobWidth);
        生成HTML(did, startTime, endTime);
        //查看速度(did, startTime, endTime, LocalDateTimeUtil.parse("20251024065125", "yyyyMMddHHmmss"), LocalDateTimeUtil.parse("20251024071818", "yyyyMMddHHmmss"));
    }

    @Test
    void 测试10秒间隔001() {
        // 一块地，非常多的路
        String did = "EM9101B8F5AZT0041";
        String yyyyMMdd = "20251027";
        String startTime = yyyyMMdd + "000000";
        String endTime = yyyyMMdd + "235959";
        double jobWidth = 3.5;
        生成数据文件(did, startTime, endTime);
        测试拆分数据(did, startTime, endTime, jobWidth);
        生成HTML(did, startTime, endTime);
    }

    @Test
    void 测试10秒间隔002() {
        // 五块地，非常分散，路也非常多，还有一个非常小的轮廓无法去掉
        String did = "EM9101B8F5AZT0041";
        String yyyyMMdd = "20251024";
        String startTime = yyyyMMdd + "000000";
        String endTime = yyyyMMdd + "235959";
        double jobWidth = 3.5;
        生成数据文件(did, startTime, endTime);
        测试拆分数据(did, startTime, endTime, jobWidth);
        生成HTML(did, startTime, endTime);
    }

    @Test
    void 测试10秒间隔003() {
        // 形成不了作业地块，不是正常的作业轨迹
        String did = "EM9101B8F5AZT0041";
        String yyyyMMdd = "20250903";
        String startTime = yyyyMMdd + "000000";
        String endTime = yyyyMMdd + "235959";
        double jobWidth = 3.5;
        生成数据文件(did, startTime, endTime);
        测试拆分数据(did, startTime, endTime, jobWidth);
        生成HTML(did, startTime, endTime);
    }

    @Test
    void 测试10秒间隔004() {
        // 应该是四块地，但是有一条轨迹也被识别为了小地块，稍微还带了点路
        String did = "EM9101B8F5AZT0041";
        String yyyyMMdd = "20251025";
        String startTime = yyyyMMdd + "000000";
        String endTime = yyyyMMdd + "235959";
        double jobWidth = 3.5;
        生成数据文件(did, startTime, endTime);
        测试拆分数据(did, startTime, endTime, jobWidth);
        生成HTML(did, startTime, endTime);
    }

    @Test
    void 测试10秒间隔005() {
        // 这个粘连了很多路，需要更多优化才行
        String did = "EM9101B8F5AZT0041";
        String yyyyMMdd = "20251026";
        String startTime = yyyyMMdd + "000000";
        String endTime = yyyyMMdd + "235959";
        double jobWidth = 3.5;
        生成数据文件(did, startTime, endTime);
        测试拆分数据(did, startTime, endTime, jobWidth);
        生成HTML(did, startTime, endTime);
    }

    @Test
    void 测试10秒间隔006() {
        String did = "EM9101B8F5AZT0041";
        String yyyyMMdd = "20251023";
        String startTime = yyyyMMdd + "000000";
        String endTime = yyyyMMdd + "235959";
        double jobWidth = 3.5;
        生成数据文件(did, startTime, endTime);
        测试拆分数据(did, startTime, endTime, jobWidth);
        生成HTML(did, startTime, endTime);
    }

    @Test
    void 测试10秒间隔007() {
        String did = "EM9101B8F4AZR0296";
        String yyyyMMdd = "20251101";
        String startTime = yyyyMMdd + "000000";
        String endTime = yyyyMMdd + "235959";
        double jobWidth = 3.5;
        生成数据文件(did, startTime, endTime);
        测试拆分数据(did, startTime, endTime, jobWidth);
        生成HTML(did, startTime, endTime);
    }


}
