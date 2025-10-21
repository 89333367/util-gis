package sunyu.util.test;

import cn.hutool.core.date.LocalDateTimeUtil;
import cn.hutool.core.io.resource.ResourceUtil;
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
import sunyu.util.GisUtil;
import sunyu.util.TDengineUtil;
import sunyu.util.pojo.TrackPoint;

import javax.sql.DataSource;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Date;
import java.util.List;
import java.util.Map;

public class TestGisUtil {
    Log log = LogFactory.get();
    GisUtil gisUtil = GisUtil.builder().build();
    ProtocolSdk protocolSdk = new ProtocolSdk("http://192.168.11.8/config.xml");

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

    private TDengineUtil getTdengineUtil() {
        return TDengineUtil.builder().dataSource(getTdengineDatasource()).build();
    }

    @Test
    void t001() throws SQLException {
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
                    log.debug("{} {} {} {} {} {} {}", did, jobArea, jobStartTime, jobEndTime, jobWidth, insertTime, updateTime);

                    String tdSql = StrUtil.format("select protocol from frequent.d_p where did='{}' and `3014`>='{}' and `3014`<='{}' and protocol match '(,2601:0,)'", did, jobStartTime, jobEndTime);
                    List<Map<String, Object>> rows = tDengineUtil.executeQuery(tdSql);
                    if (rows.isEmpty()) {
                        log.error("{} 异常，在 {} {} 找不到数据", did, jobStartTime, jobEndTime);
                    } else {
                        log.debug("找到 {} 条", rows.size());
                        List<TrackPoint> l = new ArrayList<>();
                        for (Map<String, Object> row : rows) {
                            Map<String, String> protocol = protocolSdk.parseProtocolString(row.get("protocol").toString());
                            l.add(new TrackPoint(Double.parseDouble(protocol.get("2602")), Double.parseDouble(protocol.get("2603")), LocalDateTimeUtil.parse(protocol.get("3014"), "yyyyMMddHHmmss")));
                        }
                        try {
                            Geometry g = gisUtil.buildOutline(l, jobWidth);
                            log.debug("轮廓创建完毕 {}", g);
                            String wkt = gisUtil.toWkt(g);
                            log.debug("WKT {}", wkt);
                            double wktMu = gisUtil.calcMu(wkt);
                            double mu = gisUtil.calcMu(g);
                            log.info("设备号：{} 作业时间：{} {} 宽幅：{} 原亩数：{} wkt亩数：{} 几何图形亩数：{}", did, jobStartTime, jobEndTime, jobWidth, jobArea
                                    , wktMu, mu);
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
        double mu = gisUtil.calcMu(ResourceUtil.readUtf8Str("datas/MULTIPOLYGON_1.txt"));
        log.info("{}", mu);
    }
}
