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
import cn.hutool.core.io.resource.ResourceUtil;
import cn.hutool.core.util.StrUtil;
import cn.hutool.db.Db;
import cn.hutool.db.DbUtil;
import cn.hutool.db.Entity;
import cn.hutool.log.Log;
import cn.hutool.log.LogFactory;
import cn.hutool.log.level.Level;
import sunyu.util.GisUtil;
import sunyu.util.TDengineUtil;
import sunyu.util.pojo.OutlinePart;
import sunyu.util.pojo.SplitRoadResult;
import sunyu.util.pojo.TrackPoint;

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
    void test_farm_work_split() throws SQLException {
        Db db = getMysqlDb();
        TDengineUtil tDengineUtil = getTdengineUtil();
        String sql = "select did, jobArea, jobStartTime, jobEndTime, jobWidth, insertTime, updateTime" +
                " from farm_work_split" +
                " where (did like 'NJ%' or did like 'EC%')" +
                "  and jobArea > 0" +
                " order by insertTime desc";
        int page = 1;
        int pageSize = 10;
        while (true) {
            int offset = (page - 1) * pageSize;
            List<Entity> query = db.query(StrUtil.format(sql + " limit {},{}", offset, pageSize));
            if (query.isEmpty()) {
                break;
            } else {
                for (Entity entity : query) {
                    String did = entity.getStr("did");
                    double jobArea = entity.getDouble("jobArea");
                    Date jobStartTime = entity.getDate("jobStartTime");
                    Date jobEndTime = entity.getDate("jobEndTime");
                    double jobWidth = entity.getDouble("jobWidth");
                    Date insertTime = entity.getDate("insertTime");
                    Date updateTime = entity.getDate("updateTime");
                    log.debug("{} {} {} {} {} {} {}", did, jobArea, jobStartTime, jobEndTime, jobWidth, insertTime,
                            updateTime);

                    String tdSql = StrUtil.format(
                            "select protocol from frequent.d_p where did='{}' and `3014`>='{}' and `3014`<='{}' and protocol match '(,2601:0,)'",
                            did, jobStartTime, jobEndTime);
                    List<Map<String, Object>> rows = tDengineUtil.executeQuery(tdSql);
                    if (rows.isEmpty()) {
                        log.error("{} 异常，在 {} {} 找不到数据", did, jobStartTime, jobEndTime);
                    } else {
                        log.debug("找到 {} 条", rows.size());
                        List<TrackPoint> l = new ArrayList<>();
                        for (Map<String, Object> row : rows) {
                            Map<String, String> protocol = protocolSdk
                                    .parseProtocolString(row.get("protocol").toString());
                            l.add(new TrackPoint(LocalDateTimeUtil.parse(protocol.get("3014"), "yyyyMMddHHmmss"),
                                    Double.parseDouble(protocol.get("2602")),
                                    Double.parseDouble(protocol.get("2603"))));
                        }
                        try {
                            SplitRoadResult res = gisUtil.splitRoad(l, jobWidth);
                            Geometry g = res.getOutline();
                            log.debug("轮廓创建完毕");
                            String wkt = gisUtil.toWkt(g);
                            log.debug("wkt获取完毕");
                            log.info("{}", wkt);
                            double wktMu = gisUtil.calcMu(wkt);
                            double mu = gisUtil.calcMu(g);
                            log.info("设备号：{} 作业时间：{} {} 宽幅：{} 原亩数：{} wkt亩数：{} 几何图形亩数：{}", did, jobStartTime, jobEndTime,
                                    jobWidth, jobArea, wktMu, mu);
                        } catch (Exception e) {
                            log.error(e);
                        }
                    }
                }
                page++;
            }
        }
    }

    @Test
    void wkt计算亩数() {
        double mu = gisUtil.calcMu(ResourceUtil.readUtf8Str("datas/POLYGON_1.txt"));
        log.info("{}", mu);
    }

    @Test
    void 读取farm_work() throws SQLException {
        Db db = getMysqlDb();
        /**
         * select id,
        did,
        jobArea,
        jobStartTime,
        jobEndTime,
        jobWidth,
        insertTime,
        updateTime
        from farm_work
        where jobStartTime >= '2025-05-01'
        and (did like 'NJ%' or did like 'EC%')
        and jobArea > 0
        order by insertTime desc
         */
        String sql = "select id, did, jobArea, jobStartTime, jobEndTime, jobWidth, insertTime, updateTime" +
                " from farm_work" +
                " where jobStartTime >= '2025-05-01'" +
                "  and (did like 'NJ%' or did like 'EC%')" +
                "  and jobArea > 0" +
                " order by insertTime desc";
        int page = 1;
        int pageSize = 10;
        while (true) {
            int offset = (page - 1) * pageSize;
            List<Entity> query = db.query(StrUtil.format(sql + " limit {},{}", offset, pageSize));
            if (query.isEmpty()) {
                break;
            } else {
                for (Entity entity : query) {
                    String did = entity.getStr("did");
                    double jobArea = entity.getDouble("jobArea");
                    Date jobStartTime = entity.getDate("jobStartTime");
                    Date jobEndTime = entity.getDate("jobEndTime");
                    double jobWidth = entity.getDouble("jobWidth");
                    Date insertTime = entity.getDate("insertTime");
                    Date updateTime = entity.getDate("updateTime");
                    log.debug("{} {} {} {} {} {} {}", did, jobArea, jobStartTime, jobEndTime, jobWidth, insertTime,
                            updateTime);
                    读取数据(did, DateUtil.format(jobStartTime, "yyyyMMdd"));
                    测试一天(did, DateUtil.format(jobStartTime, "yyyyMMdd"), jobWidth);
                }
                page++;
            }
        }
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

    void 读取数据(String did, String yyyyMMdd) {
        TDengineUtil tDengineUtil = getTdengineUtil();
        String jobStartTime = DateUtil.parse(yyyyMMdd, "yyyyMMdd").toString("yyyy-MM-dd") + " 00:00:00";
        String jobEndTime = DateUtil.parse(yyyyMMdd + "235959", "yyyyMMddHHmmss").toString("yyyy-MM-dd HH:mm:ss");
        String tdSql = StrUtil.format(
                "select protocol from frequent.d_p where did='{}' and `3014`>='{}' and `3014`<='{}' and protocol match '(,2601:0,)'",
                did, jobStartTime, jobEndTime);
        log.debug("{}", tdSql);
        List<Map<String, Object>> rows = tDengineUtil.executeQuery(tdSql);
        log.debug("{}", rows.size());
        List<String> l = new ArrayList<>();
        for (Map<String, Object> row : rows) {
            Map<String, String> protocol = protocolSdk.parseProtocolString(row.get("protocol").toString());
            if (Convert.toDouble(protocol.get("2602")) <= 0 || Convert.toDouble(protocol.get("2603")) <= 0) {
                continue;
            }
            l.add(StrUtil.format("{},{},{}", protocol.get("3014"), protocol.get("2602"), protocol.get("2603")));
        }
        FileUtil.writeUtf8Lines(l, path + StrUtil.format("/trace_{}_{}.txt", did, yyyyMMdd));
    }

    void 测试一天(String did, String yyyyMMdd, double jobWidth) {
        String fileName = path + StrUtil.format("/trace_{}_{}.txt", did, yyyyMMdd);
        List<TrackPoint> l = new ArrayList<>();
        DateTime jobStartTime = DateUtil.parse(yyyyMMdd, "yyyyMMdd");
        DateTime jobEndTime = DateUtil.parse(yyyyMMdd + "235959", "yyyyMMddHHmmss");
        for (String line : FileUtil.readUtf8Lines(fileName)) {
            // 20251013120625,113.33316443,28.08500825
            String[] split = line.split(",");
            TrackPoint trackPoint = new TrackPoint(LocalDateTimeUtil.parse(split[0], "yyyyMMddHHmmss"),
                    Double.parseDouble(split[1]), Double.parseDouble(split[2]));
            l.add(trackPoint);
        }
        try {
            SplitRoadResult res = gisUtil.splitRoad(l, jobWidth, 10);
            Geometry outline = res.getOutline();
            String outlineWkt = res.getWkt();
            List<OutlinePart> parts = res.getParts();

            String outlineFile = StrUtil.format(path + "/outline_{}_{}.txt", did, jobEndTime.toString("yyyyMMdd"));
            String partsFile = StrUtil.format(path + "/parts_{}_{}.txt", did, jobEndTime.toString("yyyyMMdd"));

            StringBuilder ob = new StringBuilder();
            ob.append("Outline type: ").append(outline.getGeometryType()).append('\n')
                    .append("Outline count: ")
                    .append(outline instanceof org.locationtech.jts.geom.MultiPolygon ? outline.getNumGeometries() : 1)
                    .append('\n')
                    .append("Outline mu: ").append(gisUtil.calcMu(outline)).append('\n')
                    .append("Outline WKT: ").append(outlineWkt).append('\n');
            FileUtil.writeUtf8String(ob.toString(), outlineFile);

            StringBuilder pb = new StringBuilder();
            pb.append("Parts size: ").append(parts == null ? 0 : parts.size()).append('\n');
            if (parts != null) {
                for (int i = 0; i < parts.size(); i++) {
                    OutlinePart p = parts.get(i);
                    pb.append("== Part ").append(i).append(" ==\n")
                            .append("type: ").append(p.getPolygon().getGeometryType()).append('\n')
                            .append("mu: ").append(p.getMu()).append('\n')
                            .append("startTime: ").append(p.getStartTime()).append('\n')
                            .append("endTime: ").append(p.getEndTime()).append('\n')
                            .append("wkt: ").append(p.getWkt()).append('\n');
                }
            }
            FileUtil.writeUtf8String(pb.toString(), partsFile);

            double wktMu = gisUtil.calcMu(outlineWkt);
            double mu = gisUtil.calcMu(outline);
            log.info("设备号：{} 作业时间：{} {} 宽幅：{} wkt亩数：{} 几何图形亩数：{}", did, jobStartTime, jobEndTime, jobWidth, wktMu, mu);
            log.info("结果文件已生成：outline={}, parts={}", outlineFile, partsFile);
        } catch (Exception e) {
            log.error(e);
        }
    }
}