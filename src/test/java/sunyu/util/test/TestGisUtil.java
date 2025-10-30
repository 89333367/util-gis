package sunyu.util.test;

import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Date;
import java.util.List;
import java.util.Map;

import javax.sql.DataSource;

import org.junit.jupiter.api.Test;
import org.locationtech.jts.geom.Geometry;

import com.zaxxer.hikari.HikariConfig;
import com.zaxxer.hikari.HikariDataSource;

import cn.hutool.core.convert.Convert;
import cn.hutool.core.date.DateTime;
import cn.hutool.core.date.DateUtil;
import cn.hutool.core.date.LocalDateTimeUtil;
import cn.hutool.core.io.FileUtil;
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
    ProtocolSdk protocolSdk = new ProtocolSdk("http://192.168.11.8/config.xml");
    String path = "D:/tmp/java道路拆分算法测试";

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
    void 测试两点距离() {
        CoordinatePoint p1 = new CoordinatePoint(116.55470301, 40.21296700);
        CoordinatePoint p2 = new CoordinatePoint(116.55473883, 40.21364248);
        double distance = gisUtil.haversine(p1, p2);
        log.info("{} 米", distance);
    }

    @Test
    void 循环() throws SQLException {
        Db db = getMysqlDb();
        List<Entity> rows = db.query(
                "select id,did,jobArea,jobStartTime,jobEndTime,jobWidth from farm_work where jobStartTime >= '2025-05-01' and jobStartTime < '2025-06-01' and (did like 'NJ%' or did like 'EC%') and jobArea > 0 order by insertTime desc limit 10");
        for (Entity row : rows) {
            String did = row.getStr("did");
            Date jobStartTime = row.getDate("jobStartTime");
            double jobWidth = row.getDouble("jobWidth");
            String yyyyMMdd = DateUtil.format(jobStartTime, "yyyyMMdd");
            读取数据(did, yyyyMMdd);
            测试一天(did, yyyyMMdd, jobWidth);
        }
    }

    @Test
    void wkt计算亩数() {
        String wkt = FileUtil.readUtf8String(path + "/wkt118.txt");
        double mu = gisUtil.calcMu(wkt);
        log.info("{}", mu);

        Geometry g = gisUtil.fromWkt(wkt);
        mu = gisUtil.calcMu(g);
        log.info("{}", mu);
    }

    @Test
    void 计算轮廓() throws Exception {
        String trackStr = FileUtil.readUtf8String(path + "/track.txt");
        List<TrackPoint> seg = new ArrayList<>();
        for (String split : trackStr.split("#")) {
            String[] ss = split.split(",");
            seg.add(new TrackPoint(LocalDateTimeUtil.parse(ss[2], "yyyyMMddHHmmss"), Convert.toDouble(ss[0]),
                    Convert.toDouble(ss[1])));
        }
        OutlinePart r = gisUtil.getOutline(seg, 5);
        FileUtil.writeUtf8String(StrUtil.format("wkt: {}\nmu: {}", r.getWkt(), r.getMu()), path + "/outline.txt");
    }

    @Test
    void 计算轮廓335() throws Exception {
        String trackStr = FileUtil.readUtf8String(path + "/track335.txt");
        List<TrackPoint> seg = new ArrayList<>();
        for (String split : trackStr.split("#")) {
            String[] ss = split.split(",");
            seg.add(new TrackPoint(LocalDateTimeUtil.parse(ss[2], "yyyyMMddHHmmss"), Convert.toDouble(ss[0]),
                    Convert.toDouble(ss[1])));
        }
        OutlinePart r = gisUtil.getOutline(seg, 5);
        FileUtil.writeUtf8String(StrUtil.format("wkt: {}\nmu: {}", r.getWkt(), r.getMu()), path + "/outline335.txt");
    }

    @Test
    void 计算轮廓118() throws Exception {
        String trackStr = FileUtil.readUtf8String(path + "/track118.txt");
        List<TrackPoint> seg = new ArrayList<>();
        for (String split : trackStr.split("#")) {
            String[] ss = split.split(",");
            seg.add(new TrackPoint(LocalDateTimeUtil.parse(ss[2], "yyyyMMddHHmmss"), Convert.toDouble(ss[0]),
                    Convert.toDouble(ss[1])));
        }
        OutlinePart r = gisUtil.getOutline(seg, 5);
        FileUtil.writeUtf8String(StrUtil.format("wkt: {}\nmu: {}", r.getWkt(), r.getMu()), path + "/outline118.txt");
    }

    @Test
    void 测试相交() throws Exception {
        String wktA = FileUtil.readUtf8String(path + "/wkt118.txt");
        String wktB = FileUtil.readUtf8String(path + "/wkt335.txt");
        WktIntersectionResult r = gisUtil.intersection(wktA, wktB);
        FileUtil.writeUtf8String(StrUtil.format("wkt: {}\nmu: {}", r.getWkt(), r.getMu()), path + "/intersection.txt");
    }

    @Test
    void t001() {
        String did = "EC73BD2509060248";
        String yyyyMMdd = "20251023";
        读取数据(did, yyyyMMdd);
        测试一天(did, yyyyMMdd, 2.5);
    }

    @Test
    void t002() {
        String did = "EC73BD2508220055";
        String yyyyMMdd = "20251013";
        读取数据(did, yyyyMMdd);
        测试一天(did, yyyyMMdd, 2.6);
    }

    @Test
    void t003() {
        String did = "EC71BT2402000001";
        String yyyyMMdd = "20251015";
        读取数据(did, yyyyMMdd);
        测试一天(did, yyyyMMdd, 3);
    }

    @Test
    void t004() {
        String did = "EC73BD2509060335";
        String yyyyMMdd = "20251015";
        读取数据(did, yyyyMMdd);
        测试一天(did, yyyyMMdd, 3);
    }

    void 读取数据(String did, String yyyyMMdd) {
        if (!FileUtil.exist(path + StrUtil.format("/{}_{}_trace.txt", did, yyyyMMdd))) {
            TDengineUtil tDengineUtil = getTdengineUtil();
            String jobStartTime = DateUtil.parse(yyyyMMdd, "yyyyMMdd").toString("yyyy-MM-dd") + " 00:00:00";
            String jobEndTime = DateUtil.parse(yyyyMMdd + "235959", "yyyyMMddHHmmss").toString("yyyy-MM-dd HH:mm:ss");
            String tdSql = StrUtil.format(
                    "select protocol from frequent.d_p where did='{}' and _rowts>='{}' and _rowts<='{}'",
                    did, jobStartTime, jobEndTime);
            log.debug("{}", tdSql);
            List<Map<String, Object>> rows = tDengineUtil.executeQuery(tdSql);
            List<String> l = new ArrayList<>();
            for (Map<String, Object> row : rows) {
                Map<String, String> protocol = protocolSdk.parseProtocolString(row.get("protocol").toString());
                if (protocol.containsKey("2601")) {//定位状态,0已定位，1未定位
                    if (!Convert.toStr(protocol.get("2601"), "1").equals("0")) {
                        continue;
                    }
                }
                double lon = Convert.toDouble(protocol.get("2602"), 0.0);//经度
                double lat = Convert.toDouble(protocol.get("2603"), 0.0);//纬度
                if (lon == 0 || lat == 0) {//0也算异常的点
                    continue;
                }
                if (Math.abs(lon) > 180 || Math.abs(lat) > 90) {//经度范围不在[-180,180],纬度范围不在[-90,90]，就是异常点
                    continue;
                }
                if (protocol.containsKey("3020")) {//终端ACC状态,0关闭，1开启
                    if (!Convert.toStr(protocol.get("3020"), "0").equals("1")) {
                        continue;
                    }
                }
                if (protocol.containsKey("4031")) {//作业标识,1作业,0非作业,2暂停
                    if (!Convert.toStr(protocol.get("4031"), "0").equals("1")) {
                        continue;
                    }
                }
                l.add(StrUtil.format("{},{},{}", protocol.get("3014"), protocol.get("2602"), protocol.get("2603")));
            }
            FileUtil.writeUtf8Lines(l, path + StrUtil.format("/{}_{}_trace.txt", did, yyyyMMdd));
        }
    }

    void 测试一天(String did, String yyyyMMdd, double jobWidth) {
        String fileName = path + StrUtil.format("/{}_{}_trace.txt", did, yyyyMMdd);
        List<TrackPoint> l = new ArrayList<>();
        DateTime jobEndTime = DateUtil.parse(yyyyMMdd + "235959", "yyyyMMddHHmmss");
        for (String line : FileUtil.readUtf8Lines(fileName)) {
            // 20251013120625,113.33316443,28.08500825
            String[] split = line.split(",");
            TrackPoint trackPoint = new TrackPoint(LocalDateTimeUtil.parse(split[0], "yyyyMMddHHmmss"),
                    Double.parseDouble(split[1]), Double.parseDouble(split[2]));
            l.add(trackPoint);
        }
        try {
            SplitRoadResult res = gisUtil.splitRoad(l, jobWidth, Integer.MAX_VALUE);
            Geometry outline = res.getOutline();
            String outlineWkt = res.getWkt();
            List<OutlinePart> parts = res.getParts();

            String outlineFile = StrUtil.format(path + "/{}_{}_outline.txt", did, jobEndTime.toString("yyyyMMdd"));
            String partsFile = StrUtil.format(path + "/{}_{}_parts.txt", did, jobEndTime.toString("yyyyMMdd"));

            StringBuilder ob = new StringBuilder();
            ob.append("Outline type: ").append(outline.getGeometryType()).append('\n')
                    .append("Parts: ")
                    .append(outline instanceof org.locationtech.jts.geom.MultiPolygon ? outline.getNumGeometries() : 1)
                    .append('\n')
                    .append("mu: ").append(gisUtil.calcMu(outline)).append('\n')
                    .append("totalWidthM: ").append(res.getTotalWidthM()).append('\n')
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
}