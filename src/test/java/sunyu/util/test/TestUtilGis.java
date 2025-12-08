package sunyu.util.test;

import cn.hutool.core.convert.Convert;
import cn.hutool.core.date.DateUtil;
import cn.hutool.core.date.LocalDateTimeUtil;
import cn.hutool.core.io.FileUtil;
import cn.hutool.core.io.resource.ResourceUtil;
import cn.hutool.core.util.StrUtil;
import cn.hutool.db.Db;
import cn.hutool.db.DbUtil;
import cn.hutool.log.Log;
import cn.hutool.log.LogFactory;
import cn.hutool.log.level.Level;
import com.zaxxer.hikari.HikariConfig;
import com.zaxxer.hikari.HikariDataSource;
import org.junit.jupiter.api.Test;
import sunyu.util.GisUtil;
import sunyu.util.TDengineUtil;
import sunyu.util.pojo.GaussPoint;
import sunyu.util.pojo.Part;
import sunyu.util.pojo.SplitResult;
import sunyu.util.pojo.Wgs84Point;

import javax.sql.DataSource;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

public class TestUtilGis {
    Log log = LogFactory.get();
    GisUtil gisUtil = GisUtil.builder().build();
    ProtocolSdk protocolSdk = new ProtocolSdk(FileUtil.getInputStream(FileUtil.file("d:/tmp/config.xml")));
    String path = "D:/GitLab/util-gis/testFiles";

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

    void 生成数据文件(String did, String startTime, String endTime) {
        if (!FileUtil.exist(path + StrUtil.format("/{}_{}_{}_trace.txt", did, startTime, endTime))) {
            TDengineUtil tDengineUtil = getTdengineUtil();
            String jobStartTime = DateUtil.parse(startTime, "yyyyMMddHHmmss").toString("yyyy-MM-dd HH:mm:ss");
            String jobEndTime = DateUtil.parse(endTime, "yyyyMMddHHmmss").toString("yyyy-MM-dd HH:mm:ss");
            String tdSql = StrUtil.format("select protocol from frequent.d_p where did='{}' and _rowts>='{}' and _rowts<='{}'", did, jobStartTime, jobEndTime);
            log.debug("{}", tdSql);
            List<Map<String, Object>> rows = tDengineUtil.executeQuery(tdSql);
            List<String> traceList = new ArrayList<>();
            List<String> protocolList = new ArrayList<>();
            for (Map<String, Object> row : rows) {
                Map<String, String> protocol = protocolSdk.parseProtocolString(row.get("protocol").toString());
                if (!protocol.containsKey("3014") || !protocol.containsKey("2602") || !protocol.containsKey("2603")) {
                    continue;
                }
                protocolList.add(row.get("protocol").toString());//只要有定位时间和经纬度信息就添加到协议列表
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
            Map<String, String> protocol = protocolSdk.parseProtocolString(protocolStr);
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
        partsInfo.add(StrUtil.format("WKT: {}", splitResult.getWkt()));
        partsInfo.add(StrUtil.format("作业总面积（亩）: {}\n", splitResult.getMu()));
        for (Part part : splitResult.getParts()) {
            List<String> partInfo = new ArrayList<>();
            partInfo.add(StrUtil.format("WKT: {}", part.getWkt()));
            partInfo.add(StrUtil.format("作业面积（亩）: {}", part.getMu()));
            partInfo.add(StrUtil.format("作业时间范围: {} - {}", part.getStartTime(), part.getEndTime()));
            partsInfo.add(StrUtil.join("\n", partInfo) + "\n");
        }
        FileUtil.writeUtf8Lines(partsInfo, partsFile);

        fileName = path + StrUtil.format("/{}_{}_{}_gauss.txt", did, startTime, endTime);
        if (!FileUtil.exist(fileName)) {
            List<GaussPoint> gaussPointList = gisUtil.toGaussPointList(l);
            List<String> gaussXyList = new ArrayList<>();
            for (GaussPoint gaussPoint : gaussPointList) {
                gaussXyList.add(StrUtil.format("{},{}", gaussPoint.getGaussX(), gaussPoint.getGaussY()));
            }
            FileUtil.writeUtf8Lines(gaussXyList, fileName);
        }
    }

    void 生成HTML(String did, String startTime, String endTime) {
        // 利用 showGeometryTemplate.html 当做模版，将trace输出到轨迹的TAB中
        String html = ResourceUtil.readUtf8Str("showGeometryTemplate.html");

        String fileName = path + StrUtil.format("/{}_{}_{}_trace.txt", did, startTime, endTime);
        if (!FileUtil.exist(fileName)) {
            return;
        }
        String trace = FileUtil.readUtf8String(fileName);
        html = StrUtil.replace(html, "${trace}", trace);

        String outline = FileUtil.readUtf8Lines(path + StrUtil.format("/{}_{}_{}_parts.txt", did, startTime, endTime)).get(1);
        html = StrUtil.replace(html, "${outline}", outline.replace("WKT: ", ""));

        FileUtil.writeUtf8String(html, path + StrUtil.format("/{}_{}_{}.html", did, startTime, endTime));
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
            Map<String, String> protocol = protocolSdk.parseProtocolString(protocolStr);
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
    void 测试1秒间隔001() {
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
        String did = "EC71BT2406060220";
        String startTime = "20251102130028";
        String endTime = "20251102153804";
        double jobWidth = 1.75;
        生成数据文件(did, startTime, endTime);
        测试拆分数据(did, startTime, endTime, jobWidth);
        生成HTML(did, startTime, endTime);
    }

    @Test
    void 测试10秒间隔001() {
        String did = "EM9101B8F5AZT0041";
        String yyyyMMdd = "20251027";
        String startTime = yyyyMMdd + "000000";
        String endTime = yyyyMMdd + "235959";
        double jobWidth = 2.5;
        生成数据文件(did, startTime, endTime);
        测试拆分数据(did, startTime, endTime, jobWidth);
        生成HTML(did, startTime, endTime);
    }
}
