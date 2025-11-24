package sunyu.util.test;

import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import javax.sql.DataSource;

import org.junit.jupiter.api.Test;
import org.locationtech.jts.geom.Geometry;

import com.zaxxer.hikari.HikariConfig;
import com.zaxxer.hikari.HikariDataSource;

import cn.hutool.core.collection.CollUtil;
import cn.hutool.core.convert.Convert;
import cn.hutool.core.date.DateField;
import cn.hutool.core.date.DateTime;
import cn.hutool.core.date.DateUtil;
import cn.hutool.core.date.LocalDateTimeUtil;
import cn.hutool.core.io.FileUtil;
import cn.hutool.core.io.resource.ResourceUtil;
import cn.hutool.core.text.csv.CsvUtil;
import cn.hutool.core.util.StrUtil;
import cn.hutool.db.Db;
import cn.hutool.db.DbUtil;
import cn.hutool.db.Entity;
import cn.hutool.log.Log;
import cn.hutool.log.LogFactory;
import cn.hutool.log.level.Level;
import sunyu.util.GisUtil;
import sunyu.util.TDengineUtil;
import sunyu.util.pojo.CoordinatePoint;
import sunyu.util.pojo.OutlinePart;
import sunyu.util.pojo.SplitRoadResult;
import sunyu.util.pojo.TrackPoint;
import sunyu.util.pojo.WktIntersectionResult;

public class TestGisUtil {
    Log log = LogFactory.get();
    GisUtil gisUtil = GisUtil.builder().build();
    ProtocolSdk protocolSdk = new ProtocolSdk(FileUtil.getInputStream(FileUtil.file("d:/tmp/config.xml")));
    String path = "D:/GitLab/util-gis/testFiles";

    private DataSource getMySqlDatasource() {
        HikariConfig config = new HikariConfig();
        config.setJdbcUrl(
                "jdbc:mysql://172.16.1.59:3306/farm?useUnicode=true&characterEncoding=UTF-8&serverTimezone=UTC&zeroDateTimeBehavior=convertToNull&useInformationSchema=true&useSSL=false&allowMultiQueries=true");
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

    private TDengineUtil getTdengineUtil() {
        return TDengineUtil.builder().dataSource(getTdengineDatasource()).build();
    }

    @Test
    void 测试点是否在圆中() {
        CoordinatePoint p = new CoordinatePoint(100.401807, 23.443696);
        CoordinatePoint center = new CoordinatePoint(100.27786, 23.60424);
        double radius = 1000.0;//米
        boolean isIn = gisUtil.inCircle(p, center, radius);
        log.info("点是否在圆中: {}", isIn);
    }

    @Test
    void 测试点是否在矩形中() {
        CoordinatePoint p = new CoordinatePoint(116.55470301, 40.21296700);
        CoordinatePoint p1 = new CoordinatePoint(116.55560000, 40.21296700); // 向东偏移约100米
        CoordinatePoint p2 = new CoordinatePoint(116.55560000, 40.21364248); // 向北偏移约100米
        boolean isIn = gisUtil.inRectangle(p, p1, p2);
        log.info("点是否在矩形中: {}", isIn);
    }

    @Test
    void 测试点是否在多边形中() {
        Geometry geom = gisUtil.toWgs84Geometry(
                "POLYGON((116.55470301 40.21296700, 116.55560000 40.21296700, 116.55560000 40.21364248, 116.55470301 40.21364248, 116.55470301 40.21296700))");

        // 测试多边形内部的点
        CoordinatePoint p1 = new CoordinatePoint(116.55515000, 40.21330000);
        boolean isIn1 = gisUtil.inGeometry(p1, geom);
        log.info("内部点[116.55515000, 40.21330000]是否在多边形中: {}", isIn1);

        // 测试多边形顶点（边界点）
        CoordinatePoint p2 = new CoordinatePoint(116.55560000, 40.21364248);
        boolean isIn2 = gisUtil.inGeometry(p2, geom);
        log.info("顶点[116.55560000, 40.21364248]是否在多边形中: {}", isIn2);

        // 测试多边形外部的点
        CoordinatePoint p3 = new CoordinatePoint(116.55600000, 40.21400000);
        boolean isIn3 = gisUtil.inGeometry(p3, geom);
        log.info("外部点[116.55600000, 40.21400000]是否在多边形中: {}", isIn3);
    }

    @Test
    void 测试两点距离() {
        /* CoordinatePoint p1 = new CoordinatePoint(116.55470301, 40.21296700);
        CoordinatePoint p2 = new CoordinatePoint(116.55473883, 40.21364248); */
        CoordinatePoint p1 = new CoordinatePoint(100.401807, 23.443696);
        CoordinatePoint p2 = new CoordinatePoint(100.27786, 23.60424);
        double distance = gisUtil.haversine(p1, p2);
        log.info("{} 米", distance);
    }

    @Test
    void 测试镂空作业轮廓1() throws Exception {
        String did = "EC71BT2406060220";
        String startTime = "20251102154200";
        String endTime = "20251102172202";
        double jobWidth = 1.75;
        读取一段轨迹数据(did, startTime, endTime);
        测试一段拆分数据(did, startTime, endTime, jobWidth);
        输出一段HTML(did, startTime, endTime);
    }

    @Test
    void 测试镂空作业轮廓2() throws Exception {
        String did = "EC71BT2406060220";
        String startTime = "20251102130028";
        String endTime = "20251102153804";
        double jobWidth = 1.75;
        读取一段轨迹数据(did, startTime, endTime);
        测试一段拆分数据(did, startTime, endTime, jobWidth);
        输出一段HTML(did, startTime, endTime);
    }

    @Test
    void 测试镂空作业轮廓3() throws Exception {
        String did = "EC73BD2509060398";
        String startTime = "20251103130852";
        String endTime = "20251103151309";
        double jobWidth = 1.0;
        读取一段轨迹数据(did, startTime, endTime);
        测试一段拆分数据(did, startTime, endTime, jobWidth);
        输出一段HTML(did, startTime, endTime);
    }

    @Test
    void 测试镂空作业轮廓4() throws Exception {
        String did = "EC71BT2406060220";
        String startTime = "20251103102528";
        String endTime = "20251103150242";
        double jobWidth = 1.75;
        读取一段轨迹数据(did, startTime, endTime);
        测试一段拆分数据(did, startTime, endTime, jobWidth);
        输出一段HTML(did, startTime, endTime);
    }

    @Test
    void 测试0018() throws Exception {
        String did = "EC73BD2506050018";
        String startTime = "20251104090717";
        String endTime = "20251104092257";
        double jobWidth = 2.8;
        读取一段轨迹数据(did, startTime, endTime);
        测试一段拆分数据(did, startTime, endTime, jobWidth);
        输出一段HTML(did, startTime, endTime);
    }

    @Test
    void 测试1335() throws Exception {
        String did = "EC73BD2509061335";
        String startTime = "20251104100606";
        String endTime = "20251104101419";
        double jobWidth = 2.8;
        读取一段轨迹数据(did, startTime, endTime);
        测试一段拆分数据(did, startTime, endTime, jobWidth);
        输出一段HTML(did, startTime, endTime);
    }

    @Test
    void 计算重复亩数0018_1335() throws Exception {
        String wkt1 = FileUtil.readUtf8Lines(path +
                "/EC73BD2506050018_20251104090717_20251104092257_outline.txt")
                .get(6).replace("WKT: ", "");
        String wkt2 = FileUtil.readUtf8Lines(path +
                "/EC73BD2509061335_20251104100606_20251104101419_outline.txt")
                .get(6).replace("WKT: ", "");
        WktIntersectionResult r = gisUtil.intersection(wkt1, wkt2);
        FileUtil.writeUtf8String(r.getWkt(), path +
                "/EC73BD2509061335_20251104100606_20251104101419_intersection.txt");
        FileUtil.writeUtf8String(r.getWkt(), path +
                "/EC73BD2509061335_20251104100606_20251104101419_intersection.txt");
        log.info("相交面积：{} 亩", r.getMu());
    }

    @Test
    void t001() {
        String did = "EC73BD2509060248";
        String yyyyMMdd = "20251023";
        读取一天轨迹数据(did, yyyyMMdd);
        测试一天拆分数据(did, yyyyMMdd, 2.5);
        输出一天HTML(did, yyyyMMdd);
    }

    @Test
    void t002() {
        String did = "EC73BD2508220055";
        String yyyyMMdd = "20251013";
        读取一天轨迹数据(did, yyyyMMdd);
        测试一天拆分数据(did, yyyyMMdd, 2.6);
        输出一天HTML(did, yyyyMMdd);
    }

    @Test
    void t003() {
        String did = "NJ4GBQSAX0000687";
        String yyyyMMdd = "20250503";
        读取一天轨迹数据(did, yyyyMMdd);
        测试一天拆分数据(did, yyyyMMdd, 2.7);
        输出一天HTML(did, yyyyMMdd);
    }

    @Test
    void t004() {
        String did = "EC73BD2509061335";
        String startTime = "20251029101603";
        String endTime = "20251029102521";
        double widthM = 2;
        读取一段轨迹数据(did, startTime, endTime);
        测试一段拆分数据(did, startTime, endTime, widthM);
        输出一段HTML(did, startTime, endTime);
    }

    @Test
    void t005() {
        String did = "NJ4GNBSAX0000693";
        String startTime = "20250505080044";
        String endTime = "20250505173658";
        double widthM = 2.5;
        读取一段轨迹数据(did, startTime, endTime);
        测试一段拆分数据(did, startTime, endTime, widthM);
        输出一段HTML(did, startTime, endTime);
    }

    @Test
    void t006() {
        String did = "NJ4GNBSAX0000693";
        String startTime = "20250507073041";
        String endTime = "20250507162457";
        double widthM = 2.5;
        读取一段轨迹数据(did, startTime, endTime);
        测试一段拆分数据(did, startTime, endTime, widthM);
        输出一段HTML(did, startTime, endTime);
    }

    @Test
    void t007() {
        String did = "NJ4GNBSAX0000693";
        String startTime = "20250509092721";
        String endTime = "20250509111620";
        double widthM = 2.5;
        读取一段轨迹数据(did, startTime, endTime);
        测试一段拆分数据(did, startTime, endTime, widthM);
        输出一段HTML(did, startTime, endTime);
    }

    @Test
    void t008() {
        String did = "EM9101B8F5AZT0041";
        String yyyyMMdd = "20251027";
        读取一天轨迹数据(did, yyyyMMdd);
        测试一天拆分数据(did, yyyyMMdd, 2.5);
        输出一天HTML(did, yyyyMMdd);
    }

    @Test
    void t009() {
        String did = "EM9101B8F5AZT0041";
        String yyyyMMdd = "20251024";
        读取一天轨迹数据(did, yyyyMMdd);
        测试一天拆分数据(did, yyyyMMdd, 2.5);
        输出一天HTML(did, yyyyMMdd);
    }

    @Test
    void t010() {
        String did = "EM9101B8F5AZT0041";
        String yyyyMMdd = "20250903";
        读取一天轨迹数据(did, yyyyMMdd);
        测试一天拆分数据(did, yyyyMMdd, 2.5);
        输出一天HTML(did, yyyyMMdd);
    }

    @Test
    void t011() {
        String did = "EM9101B8F5AZT0041";
        String yyyyMMdd = "20251025";
        读取一天轨迹数据(did, yyyyMMdd);
        测试一天拆分数据(did, yyyyMMdd, 2.5);
        输出一天HTML(did, yyyyMMdd);
    }

    @Test
    void t012() {
        String did = "EM9101B8F5AZT0041";
        String yyyyMMdd = "20251026";
        读取一天轨迹数据(did, yyyyMMdd);
        测试一天拆分数据(did, yyyyMMdd, 2.5);
        输出一天HTML(did, yyyyMMdd);
    }

    @Test
    void t013() {
        String did = "EM9101B8F5AZT0041";
        String yyyyMMdd = "20251023";
        读取一天轨迹数据(did, yyyyMMdd);
        测试一天拆分数据(did, yyyyMMdd, 2.5);
        输出一天HTML(did, yyyyMMdd);
    }

    @Test
    void t014() {
        String did = "EM9101B8F4AZR0296";
        String yyyyMMdd = "20251101";
        读取一天轨迹数据(did, yyyyMMdd);
        测试一天拆分数据(did, yyyyMMdd, 2.5);
        输出一天HTML(did, yyyyMMdd);
    }

    @Test
    void t015() {
        String did = "EC71BT2404140062";
        String yyyyMMdd = "20240426";
        读取一天轨迹数据(did, yyyyMMdd);
        测试一天拆分数据(did, yyyyMMdd, 2.6);
        输出一天HTML(did, yyyyMMdd);
    }

    @Test
    void t016() {
        String did = "EC71BT2404140062";
        String yyyyMMdd = "20240427";
        读取一天轨迹数据(did, yyyyMMdd);
        测试一天拆分数据(did, yyyyMMdd, 2.1);
        输出一天HTML(did, yyyyMMdd);
    }

    @Test
    void t017() {
        String did = "EC71BT2404140062";
        String yyyyMMdd = "20240428";
        读取一天轨迹数据(did, yyyyMMdd);
        测试一天拆分数据(did, yyyyMMdd, 3.9);
        输出一天HTML(did, yyyyMMdd);
    }

    @Test
    void t018() {
        String did = "EC71BT2404140062";
        String yyyyMMdd = "20240429";
        读取一天轨迹数据(did, yyyyMMdd);
        测试一天拆分数据(did, yyyyMMdd, 3.9);
        输出一天HTML(did, yyyyMMdd);
    }

    @Test
    void 循环金大丰1() {
        List<String> l = new ArrayList<>();
        String did = "EM9101B8F5AZT0041";
        String startTime = "20250903000000";
        String endTime = "20251028000000";
        List<DateTime> dates = DateUtil.rangeToList(DateUtil.parse(startTime, "yyyyMMddHHmmss"),
                DateUtil.parse(endTime, "yyyyMMddHHmmss"), DateField.DAY_OF_YEAR);
        for (DateTime date : dates) {
            String yyyyMMdd = DateUtil.format(date, "yyyyMMdd");
            log.info("{}", yyyyMMdd);
            读取一天轨迹数据(did, yyyyMMdd);
            l.add(测试一天拆分数据(did, yyyyMMdd, 2.5));
            输出一天HTML(did, yyyyMMdd);
        }
        FileUtil.writeUtf8Lines(l, path + "/" + did + "_mu.txt");
    }

    @Test
    void 循环金大丰2() {
        List<String> l = new ArrayList<>();
        String did = "EM9101B8F4AZR0296";
        String startTime = "20251028000000";
        String endTime = "20251105000000";
        List<DateTime> dates = DateUtil.rangeToList(DateUtil.parse(startTime, "yyyyMMddHHmmss"),
                DateUtil.parse(endTime, "yyyyMMddHHmmss"), DateField.DAY_OF_YEAR);
        for (DateTime date : dates) {
            String yyyyMMdd = DateUtil.format(date, "yyyyMMdd");
            log.info("{}", yyyyMMdd);
            读取一天轨迹数据(did, yyyyMMdd);
            l.add(测试一天拆分数据(did, yyyyMMdd, 2.5));
            输出一天HTML(did, yyyyMMdd);
        }
        FileUtil.writeUtf8Lines(l, path + "/" + did + "_mu.txt");
    }

    @Test
    void 导出csv() {
        List<String> lines = FileUtil
                .readUtf8Lines("D:\\GitLab\\util-gis\\testFiles\\EM9101B8F5AZT0041_20251026_trace.txt");
        List<List<String>> rows = new ArrayList<>();
        rows.add(Arrays.asList("lon", "lat", "timestamp"));
        rows.addAll(lines.stream().map(line -> StrUtil.split(line, ",")).collect(Collectors.toList()));
        CsvUtil.getWriter("d:/tmp/EM9101B8F5AZT0041_20251026_trace.csv", null).write(rows);
    }

    @Test
    void 轨迹围栏测试() {
        String did = "BAD3322508600098";
        String yyyyMMdd = "20251124";
        读取一天轨迹数据(did, yyyyMMdd);
        测试一天拆分数据(did, yyyyMMdd, 2.5);
        输出一天HTML(did, yyyyMMdd);
    }

    @Test
    void 循环一天() throws SQLException {
        Db db = getMysqlDb();
        List<Entity> rows = db.query(
                "select id,did,jobArea,jobStartTime,jobEndTime,jobWidth from farm_work where jobStartTime >= '2025-05-01' and jobStartTime < '2025-06-01' and (did like 'NJ%' or did like 'EC%') and jobArea > 0 and fixFlag = 1 and effectiveJobArea >0 order by insertTime desc limit 10");
        for (Entity row : rows) {
            String did = row.getStr("did");
            Date jobStartTime = row.getDate("jobStartTime");
            double jobWidth = row.getDouble("jobWidth");
            String yyyyMMdd = DateUtil.format(jobStartTime, "yyyyMMdd");
            读取一天轨迹数据(did, yyyyMMdd);
            测试一天拆分数据(did, yyyyMMdd, jobWidth);
            输出一天HTML(did, yyyyMMdd);
        }
    }

    @Test
    void 循环farm_work() throws SQLException {
        Db db = getMysqlDb();
        List<Entity> rows = db.query(
                "select id,did,jobArea,jobStartTime,jobEndTime,jobWidth from farm_work where jobStartTime >= '2025-05-01' and jobStartTime < '2025-06-01' and (did like 'NJ%' or did like 'EC%') and jobArea > 0 and fixFlag = 1 and effectiveJobArea >0 order by insertTime desc limit 10");
        for (Entity row : rows) {
            String did = row.getStr("did");
            Date jobStartTime = row.getDate("jobStartTime");
            Date jobEndTime = row.getDate("jobEndTime");
            double jobWidth = row.getDouble("jobWidth");
            String startTime = DateUtil.format(jobStartTime, "yyyyMMddHHmmss");
            String endTime = DateUtil.format(jobEndTime, "yyyyMMddHHmmss");
            读取一段轨迹数据(did, startTime, endTime);
            测试一段拆分数据(did, startTime, endTime, jobWidth);
            输出一段HTML(did, startTime, endTime);
        }
    }

    void 读取一段轨迹数据(String did, String startTime, String endTime) {
        if (!FileUtil.exist(path + StrUtil.format("/{}_{}_{}_trace.txt", did, startTime, endTime))) {
            TDengineUtil tDengineUtil = getTdengineUtil();
            String jobStartTime = DateUtil.parse(startTime, "yyyyMMddHHmmss").toString("yyyy-MM-dd HH:mm:ss");
            String jobEndTime = DateUtil.parse(endTime, "yyyyMMddHHmmss").toString("yyyy-MM-dd HH:mm:ss");
            String tdSql = StrUtil.format(
                    "select protocol from frequent.d_p where did='{}' and _rowts>='{}' and _rowts<='{}'",
                    did, jobStartTime, jobEndTime);
            log.debug("{}", tdSql);
            List<Map<String, Object>> rows = tDengineUtil.executeQuery(tdSql);
            List<String> l = new ArrayList<>();
            for (Map<String, Object> row : rows) {
                Map<String, String> protocol = protocolSdk.parseProtocolString(row.get("protocol").toString());
                if (!protocol.containsKey("2601") || !protocol.containsKey("2602") || !protocol.containsKey("2603")
                        || !protocol.containsKey("3014")) {
                    continue;
                }
                if (!protocol.get("2601").equals("0")) {// 定位状态,0已定位，1未定位
                    continue;
                }
                if (protocol.containsKey("3020")) {// 终端ACC状态,0关闭，1开启
                    if (!Convert.toStr(protocol.get("3020"), "0").equals("1")) {
                        continue;
                    }
                }
                if (protocol.containsKey("4031")) {// 作业标识,1作业,0非作业,2暂停
                    if (!Convert.toStr(protocol.get("4031"), "0").equals("1")) {
                        continue;
                    }
                }
                // {定位时间yyyyMMddHHmmss},{经度},{纬度}
                l.add(StrUtil.format("{},{},{}", protocol.get("3014"), protocol.get("2602"), protocol.get("2603")));
            }
            FileUtil.writeUtf8Lines(l, path + StrUtil.format("/{}_{}_{}_trace.txt", did, startTime, endTime));
        }
    }

    void 读取一天轨迹数据(String did, String yyyyMMdd) {
        if (!FileUtil.exist(path + StrUtil.format("/{}_{}_trace.txt", did, yyyyMMdd))) {
            TDengineUtil tDengineUtil = getTdengineUtil();
            String jobStartTime = DateUtil.parse(yyyyMMdd, "yyyyMMdd").toString("yyyy-MM-dd") + " 00:00:00";
            String jobEndTime = DateUtil.parse(yyyyMMdd + "235959", "yyyyMMddHHmmss").toString("yyyy-MM-dd HH:mm:ss");
            String tdSql = StrUtil.format(
                    "select protocol from frequent.d_p where did='{}' and _rowts>='{}' and _rowts<='{}'",
                    did, jobStartTime, jobEndTime);
            log.debug("{}", tdSql);
            List<Map<String, Object>> rows = tDengineUtil.executeQuery(tdSql);
            log.info("读取到 {} 条记录", rows.size());
            List<String> l = new ArrayList<>();
            for (Map<String, Object> row : rows) {
                Map<String, String> protocol = protocolSdk.parseProtocolString(row.get("protocol").toString());
                if (!protocol.containsKey("2601") || !protocol.containsKey("2602") || !protocol.containsKey("2603")
                        || !protocol.containsKey("2204") || !protocol.containsKey("3012")
                        || !protocol.containsKey("3014")) {
                    continue;
                }
                if (!protocol.get("2601").equals("0")) {// 定位状态,0已定位，1未定位
                    continue;
                }
                if (protocol.containsKey("3020")) {// 终端ACC状态,0关闭，1开启
                    if (!Convert.toStr(protocol.get("3020"), "0").equals("1")) {
                        continue;
                    }
                }
                if (protocol.containsKey("4031")) {// 作业标识,1作业,0非作业,2暂停
                    if (!Convert.toStr(protocol.get("4031"), "0").equals("1")) {
                        continue;
                    }
                }
                // {定位时间yyyyMMddHHmmss},{经度},{纬度}
                l.add(StrUtil.format("{},{},{}", protocol.get("3014"), protocol.get("2602"), protocol.get("2603")));
            }
            if (CollUtil.isNotEmpty(l)) {
                FileUtil.writeUtf8Lines(l, path + StrUtil.format("/{}_{}_trace.txt", did, yyyyMMdd));
            }
        }
    }

    String 测试一天拆分数据(String did, String yyyyMMdd, double jobWidth) {
        String fileName = path + StrUtil.format("/{}_{}_trace.txt", did, yyyyMMdd);
        if (!FileUtil.exist(fileName)) {
            return null;
        }
        List<TrackPoint> l = new ArrayList<>();
        DateTime jobEndTime = DateUtil.parse(yyyyMMdd + "235959", "yyyyMMddHHmmss");
        for (String line : FileUtil.readUtf8Lines(fileName)) {
            String[] split = line.split(",");
            TrackPoint trackPoint = new TrackPoint();
            trackPoint.setTime(LocalDateTimeUtil.parse(split[0], "yyyyMMddHHmmss"));
            trackPoint.setLon(Double.parseDouble(split[1]));
            trackPoint.setLat(Double.parseDouble(split[2]));
            l.add(trackPoint);
        }
        try {
            SplitRoadResult res = gisUtil.splitRoad(l, jobWidth);
            Geometry outline = res.getOutline();
            String outlineWkt = res.getWkt();
            List<OutlinePart> parts = res.getParts();

            String outlineFile = StrUtil.format(path + "/{}_{}_outline.txt", did, jobEndTime.toString("yyyyMMdd"));
            String partsFile = StrUtil.format(path + "/{}_{}_parts.txt", did, jobEndTime.toString("yyyyMMdd"));

            StringBuilder ob = new StringBuilder();
            ob.append("type: ").append(outline.getGeometryType()).append('\n')
                    .append("Parts: ")
                    .append(outline instanceof org.locationtech.jts.geom.MultiPolygon ? outline.getNumGeometries() : 1)
                    .append('\n')
                    .append("mu: ").append(res.getMu()).append('\n')
                    .append("totalWidthM: ").append(res.getTotalWidthM()).append('\n')
                    .append("startTime: ").append(res.getStartTime()).append('\n')
                    .append("endTime: ").append(res.getEndTime()).append('\n')
                    .append("WKT: ").append(outlineWkt).append('\n');
            FileUtil.writeUtf8String(ob.toString(), outlineFile);

            StringBuilder pb = new StringBuilder();
            pb.append("Parts: ").append(parts == null ? 0 : parts.size()).append('\n');
            if (parts != null) {
                for (int i = 0; i < parts.size(); i++) {
                    OutlinePart p = parts.get(i);
                    pb.append("== Part ").append(i).append(" ==\n")
                            .append("type: ").append(p.getOutline().getGeometryType()).append('\n')
                            .append("mu: ").append(p.getMu()).append('\n')
                            .append("totalWidthM: ").append(p.getTotalWidthM()).append('\n')
                            .append("startTime: ").append(p.getStartTime()).append('\n')
                            .append("endTime: ").append(p.getEndTime()).append('\n')
                            .append("wkt: ").append(p.getWkt()).append('\n')
                            .append("trackStr: ").append(p.getTrackStr()).append('\n');
                }
            }
            FileUtil.writeUtf8String(pb.toString(), partsFile);

            log.info("结果文件已生成：outline={}, parts={}", outlineFile, partsFile);

            return StrUtil.format("{} {} {}", did, yyyyMMdd, res.getMu());
        } catch (Exception e) {
            log.error(e);
        }
        return null;
    }

    void 测试一段拆分数据(String did, String startTime, String endTime, double jobWidth) {
        String fileName = path + StrUtil.format("/{}_{}_{}_trace.txt", did, startTime, endTime);
        if (!FileUtil.exist(fileName)) {
            return;
        }
        List<TrackPoint> l = new ArrayList<>();
        for (String line : FileUtil.readUtf8Lines(fileName)) {
            // 20251013120625,113.33316443,28.08500825
            String[] split = line.split(",");
            TrackPoint trackPoint = new TrackPoint();
            trackPoint.setTime(LocalDateTimeUtil.parse(split[0], "yyyyMMddHHmmss"));
            trackPoint.setLon(Double.parseDouble(split[1]));
            trackPoint.setLat(Double.parseDouble(split[2]));
            l.add(trackPoint);
        }
        try {
            SplitRoadResult res = gisUtil.splitRoad(l, jobWidth);
            Geometry outline = res.getOutline();
            String outlineWkt = res.getWkt();
            List<OutlinePart> parts = res.getParts();

            String outlineFile = StrUtil.format(path + "/{}_{}_{}_outline.txt", did, startTime, endTime);
            String partsFile = StrUtil.format(path + "/{}_{}_{}_parts.txt", did, startTime, endTime);

            StringBuilder ob = new StringBuilder();
            ob.append("type: ").append(outline.getGeometryType()).append('\n')
                    .append("Parts: ")
                    .append(outline instanceof org.locationtech.jts.geom.MultiPolygon ? outline.getNumGeometries() : 1)
                    .append('\n')
                    .append("mu: ").append(res.getMu()).append('\n')
                    .append("totalWidthM: ").append(res.getTotalWidthM()).append('\n')
                    .append("startTime: ").append(res.getStartTime()).append('\n')
                    .append("endTime: ").append(res.getEndTime()).append('\n')
                    .append("WKT: ").append(outlineWkt).append('\n');
            FileUtil.writeUtf8String(ob.toString(), outlineFile);

            StringBuilder pb = new StringBuilder();
            pb.append("Parts: ").append(parts == null ? 0 : parts.size()).append('\n');
            if (parts != null) {
                for (int i = 0; i < parts.size(); i++) {
                    OutlinePart p = parts.get(i);
                    pb.append("== Part ").append(i).append(" ==\n")
                            .append("type: ").append(p.getOutline().getGeometryType()).append('\n')
                            .append("mu: ").append(p.getMu()).append('\n')
                            .append("totalWidthM: ").append(p.getTotalWidthM()).append('\n')
                            .append("startTime: ").append(p.getStartTime()).append('\n')
                            .append("endTime: ").append(p.getEndTime()).append('\n')
                            .append("wkt: ").append(p.getWkt()).append('\n')
                            .append("trackStr: ").append(p.getTrackStr()).append('\n');
                }
            }
            FileUtil.writeUtf8String(pb.toString(), partsFile);

            log.info("结果文件已生成：outline={}, parts={}", outlineFile, partsFile);
        } catch (Exception e) {
            log.error(e);
        }
    }

    void 输出一天HTML(String did, String yyyyMMdd) {
        // 利用 showGeometryTemplate.html 当做模版，将trace输出到轨迹的TAB中
        String html = ResourceUtil.readUtf8Str("showGeometryTemplate.html");

        String fileName = path + StrUtil.format("/{}_{}_trace.txt", did, yyyyMMdd);
        if (!FileUtil.exist(fileName)) {
            return;
        }
        String trace = FileUtil.readUtf8String(fileName);
        html = StrUtil.replace(html, "${trace}", trace);

        String outline = FileUtil.readUtf8Lines(path + StrUtil.format("/{}_{}_outline.txt", did, yyyyMMdd)).get(6);
        html = StrUtil.replace(html, "${outline}", outline.replace("WKT: ", ""));

        FileUtil.writeUtf8String(html, path + StrUtil.format("/{}_{}.html", did, yyyyMMdd));
    }

    void 输出一段HTML(String did, String startTime, String endTime) {
        // 利用 showGeometryTemplate.html 当做模版，将trace输出到轨迹的TAB中
        String html = ResourceUtil.readUtf8Str("showGeometryTemplate.html");

        String fileName = path + StrUtil.format("/{}_{}_{}_trace.txt", did, startTime, endTime);
        if (!FileUtil.exist(fileName)) {
            return;
        }
        String trace = FileUtil.readUtf8String(fileName);
        html = StrUtil.replace(html, "${trace}", trace);

        String outline = FileUtil.readUtf8Lines(path + StrUtil.format("/{}_{}_{}_outline.txt", did, startTime, endTime))
                .get(6);
        html = StrUtil.replace(html, "${outline}", outline.replace("WKT: ", ""));

        FileUtil.writeUtf8String(html, path + StrUtil.format("/{}_{}_{}.html", did, startTime, endTime));
    }
}