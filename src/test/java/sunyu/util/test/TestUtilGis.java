package sunyu.util.test;

import cn.hutool.core.convert.Convert;
import cn.hutool.core.date.LocalDateTimeUtil;
import cn.hutool.core.io.FileUtil;
import cn.hutool.core.io.resource.ResourceUtil;
import cn.hutool.core.util.ReUtil;
import cn.hutool.core.util.StrUtil;
import cn.hutool.log.Log;
import cn.hutool.log.LogFactory;
import org.junit.jupiter.api.Test;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.MultiPolygon;
import sunyu.util.GisUtil;
import sunyu.util.pojo.*;
import sunyu.util.test.config.MyBatis;
import sunyu.util.test.entity.DP;
import sunyu.util.test.entity.FarmWork;
import sunyu.util.test.mapper.farm.FarmMapper;
import sunyu.util.test.mapper.tdengine.TdengineMapper;

import java.time.LocalDateTime;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

public class TestUtilGis {
    Log log = LogFactory.get();
    GisUtil gisUtil = GisUtil.builder().build();
    String path = "D:/tmp/testFiles";
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

    List<DP> selectWorkPoints(String did, LocalDateTime startTime, LocalDateTime endTime) {
        TdengineMapper mapper = MyBatis.getMapper(TdengineMapper.class);
        List<DP> l = mapper.selectWorkPoints(did, startTime, endTime, false);
        if (l.isEmpty()) {
            l = mapper.selectWorkProtocol(did, startTime, endTime, false);
            for (DP dp : l) {
                Map<String, String> protocol = parseProtocolString(dp.getProtocol());
                dp.setP3014(LocalDateTimeUtil.parse(protocol.get("3014"), "yyyyMMddHHmmss"));
                dp.setP2602(Double.parseDouble(protocol.get("2602")));
                dp.setP2603(Double.parseDouble(protocol.get("2603")));
                dp.setP3020(Convert.toInt(protocol.get("3020")));
                dp.setP4031(Convert.toInt(protocol.get("4031")));
            }
        }
        return l;
    }

    void 测试拆分数据(String did, String startTime, String endTime, double jobWidth, boolean updateFarmWorkTable) {
        测试拆分数据(did, startTime, endTime, jobWidth, updateFarmWorkTable, new SplitRoadParams());
    }

    void 测试拆分数据(String did, String startTime, String endTime, double jobWidth) {
        测试拆分数据(did, startTime, endTime, jobWidth, false, new SplitRoadParams());
    }

    void 测试拆分数据(String did, String startTime, String endTime, double jobWidth, Boolean updateFarmWorkTable, SplitRoadParams splitRoadParams) {
        List<DP> dps = selectWorkPoints(did, LocalDateTimeUtil.parse(startTime, "yyyyMMddHHmmss"), LocalDateTimeUtil.parse(endTime, "yyyyMMddHHmmss"));
        List<Wgs84Point> l = new ArrayList<>();
        for (DP dp : dps) {
            Wgs84Point wgs84Point = new Wgs84Point();
            wgs84Point.setGpsTime(dp.getP3014());
            wgs84Point.setLongitude(dp.getP2602());
            wgs84Point.setLatitude(dp.getP2603());
            if (dp.getP3020() != null) {
                wgs84Point.setJobStatus(2);//先设置非作业状态
                if (dp.getP3020() == 1) {// 终端ACC状态,0关闭，1开启
                    wgs84Point.setJobStatus(1);
                }
            }
            if (splitRoadParams.getCheckWorkingStatus() && dp.getP4031() != null) {// 作业标识,1作业,0非作业,2暂停
                wgs84Point.setJobStatus(2);//先设置非作业状态
                if (dp.getP4031() == 1) {
                    wgs84Point.setJobStatus(1);//作业标识是1，认为是作业
                }
            }
            l.add(wgs84Point);
        }

        String partsFile = StrUtil.format(path + "/{}_{}_{}_parts.txt", did, startTime, endTime);
        List<String> partsInfo = new ArrayList<>();
        SplitResult splitResult = gisUtil.splitRoad(l, jobWidth, splitRoadParams);
        partsInfo.add(StrUtil.format("共有 {} 个地块", splitResult.getGaussGeometry().getNumGeometries()));
        partsInfo.add(StrUtil.format("作业总幅宽（米）: {}", splitResult.getWorkingWidth()));
        partsInfo.add(StrUtil.format("总WKT: {}", splitResult.getWkt()));
        partsInfo.add(StrUtil.format("作业总面积（亩）: {}", splitResult.getMu()));
        partsInfo.add(StrUtil.format("作业时间范围: {} - {}", splitResult.getStartTime(), splitResult.getEndTime()));
        partsInfo.add(StrUtil.format("最小有效时间间隔（秒）: {}", splitResult.getMinEffectiveInterval()));
        partsInfo.add(StrUtil.format("聚类总点数: {}", splitResult.getClusterPointCount()));
        if (splitResult.getCenterWgs84Point() != null) {
            partsInfo.add(StrUtil.format("中心点：经度{},纬度{}", splitResult.getCenterWgs84Point().getLongitude(), splitResult.getCenterWgs84Point().getLatitude()));
        }
        partsInfo.add("\n");
        int partIndex = 1;
        for (FarmPlot farmPlot : splitResult.getFarmPlots()) {
            List<String> partInfo = new ArrayList<>();
            partInfo.add(StrUtil.format("地块 {}:", partIndex++));
            if (farmPlot.getGaussGeometry() instanceof MultiPolygon) {
                partInfo.add(StrUtil.format("包含 {} 个子地块", farmPlot.getGaussGeometry().getNumGeometries()));
            }
            partInfo.add(StrUtil.format("作业总幅宽（米）：{}", farmPlot.getWorkingWidth()));
            partInfo.add(StrUtil.format("子WKT: {}", farmPlot.getWkt()));
            partInfo.add(StrUtil.format("作业面积（亩）: {}", farmPlot.getMu()));
            partInfo.add(StrUtil.format("作业时间范围: {} - {}", farmPlot.getStartTime(), farmPlot.getEndTime()));
            partInfo.add(StrUtil.format("最小有效时间间隔（秒）: {}", farmPlot.getMinEffectiveInterval()));
            partInfo.add(StrUtil.format("聚类点数: {}", farmPlot.getClusterPointCount()));
            if (farmPlot.getCenterWgs84Point() != null) {
                partInfo.add(StrUtil.format("中心点：经度{},纬度{}", farmPlot.getCenterWgs84Point().getLongitude(), farmPlot.getCenterWgs84Point().getLatitude()));
            }
            partsInfo.add(StrUtil.join("\n", partInfo) + "\n");
        }
        FileUtil.writeUtf8Lines(partsInfo, partsFile);

        String fileName = path + StrUtil.format("/{}_{}_{}_gauss.txt", did, startTime, endTime);
        l = gisUtil.filterWgs84Points(l);
        List<GaussPoint> gaussPointList = gisUtil.toGaussPointList(l);
        List<String> gaussXyList = new ArrayList<>();
        for (GaussPoint gaussPoint : gaussPointList) {
            gaussXyList.add(StrUtil.format("{},{}", gaussPoint.getGaussX(), gaussPoint.getGaussY()));
        }
        FileUtil.writeUtf8Lines(gaussXyList, fileName);

        String html = ResourceUtil.readUtf8Str("showGeometrysTemplate_leaflet.html");
        StringBuilder trace = new StringBuilder();
        for (Wgs84Point wgs84Point : l) {
            trace.append(StrUtil.format("{},{},{}\n", LocalDateTimeUtil.format(wgs84Point.getGpsTime(), "yyyyMMddHHmmss"), wgs84Point.getLongitude(), wgs84Point.getLatitude()));
        }
        html = StrUtil.replace(html, "${trace}", trace.toString());
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

        if (updateFarmWorkTable != null && updateFarmWorkTable) {
            // 直接更新farm_work表的复算亩数以及WKT
            FarmMapper mapper = MyBatis.getMapper(FarmMapper.class);
            FarmWork farmWork = new FarmWork();
            farmWork.setDid(did);
            farmWork.setJobEndTime(LocalDateTimeUtil.parse(endTime, "yyyyMMddHHmmss"));
            farmWork.setEffectiveJobArea(splitResult.getMu());
            farmWork.setWktPoly(splitResult.getWkt());
            mapper.updateFarmWork(farmWork);
        }
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
        // 中间有镂空，识别出1个地块
        String did = "EC71BT2406060220";
        String startTime = "20251102154200";
        String endTime = "20251102172202";
        double jobWidth = 1.75;
        测试拆分数据(did, startTime, endTime, jobWidth);
    }

    @Test
    void 测试1秒间隔002() {
        // 识别出3个地块
        String did = "EC71BT2406060220";
        String startTime = "20251102130028";
        String endTime = "20251102153804";
        double jobWidth = 1.75;
        测试拆分数据(did, startTime, endTime, jobWidth);
    }

    @Test
    void 测试1秒间隔003() {
        // 乱跑的，这个地块处理有些问题
        String did = "EC73BD2509060398";
        String startTime = "20251103130852";
        String endTime = "20251103151309";
        double jobWidth = 1.0;
        测试拆分数据(did, startTime, endTime, jobWidth);
    }

    @Test
    void 测试1秒间隔004() {
        // 识别出3个地块，中间有个大石头，有一块地中有一条线被多切割了一点
        String did = "EC71BT2406060220";
        String startTime = "20251103102528";
        String endTime = "20251103150242";
        double jobWidth = 1.75;
        测试拆分数据(did, startTime, endTime, jobWidth);
    }

    @Test
    void 测试1秒间隔005() {
        // 广西田测，测试地块，1.26亩
        String did = "EC73BD2506050018";
        String startTime = "20251104090717";
        String endTime = "20251104092257";
        double jobWidth = 2.8;
        测试拆分数据(did, startTime, endTime, jobWidth);
    }

    @Test
    void 测试1秒间隔006() {
        // 广西田测，测试地块，0.62亩
        String did = "EC73BD2509061335";
        String startTime = "20251104100606";
        String endTime = "20251104101419";
        double jobWidth = 2.8;
        测试拆分数据(did, startTime, endTime, jobWidth);
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
        // 识别出2个地块，有一块类似于蝴蝶形状
        String did = "EC73BD2509060248";
        String yyyyMMdd = "20251023";
        String startTime = yyyyMMdd + "000000";
        String endTime = yyyyMMdd + "235959";
        double jobWidth = 2.5;
        测试拆分数据(did, startTime, endTime, jobWidth);
    }

    @Test
    void 测试1秒间隔008() {
        // 一天轨迹，跑了多个地块，非常分散
        String did = "EC73BD2508220055";
        String yyyyMMdd = "20251013";
        String startTime = yyyyMMdd + "000000";
        String endTime = yyyyMMdd + "235959";
        double jobWidth = 2.6;
        测试拆分数据(did, startTime, endTime, jobWidth);
    }

    @Test
    void 测试1秒间隔009() {
        // 识别出两个地块
        String did = "NJ4GBQSAX0000687";
        String yyyyMMdd = "20250503";
        String startTime = yyyyMMdd + "000000";
        String endTime = yyyyMMdd + "235959";
        double jobWidth = 2.7;
        测试拆分数据(did, startTime, endTime, jobWidth);
    }

    @Test
    void 测试1秒间隔010() {
        // 不会形成作业地块，不是正常的作业轨迹
        String did = "EC73BD2509061335";
        String startTime = "20251029101603";
        String endTime = "20251029102521";
        double jobWidth = 2;
        测试拆分数据(did, startTime, endTime, jobWidth);
    }

    @Test
    void 测试1秒间隔011() {
        // 识别出了5个地块
        String did = "NJ4GNBSAX0000693";
        String startTime = "20250505080044";
        String endTime = "20250505173658";
        double jobWidth = 2.5;
        测试拆分数据(did, startTime, endTime, jobWidth);
    }

    @Test
    void 测试1秒间隔012() {
        // 识别出了3个地块，有一个非常小的块，被删掉了
        String did = "NJ4GNBSAX0000693";
        String startTime = "20250507073041";
        String endTime = "20250507162457";
        double jobWidth = 2.5;
        测试拆分数据(did, startTime, endTime, jobWidth);
    }

    @Test
    void 测试1秒间隔013() {
        // 识别出1块地
        String did = "NJ4GNBSAX0000693";
        String startTime = "20250509092721";
        String endTime = "20250509111620";
        double jobWidth = 2.5;
        //测试拆分数据(did, startTime, endTime, jobWidth);
        测试拆分数据(did, startTime, endTime, jobWidth, false, new SplitRoadParams().setRoadWidth(5.0));
    }

    @Test
    void 测试1秒间隔014() {
        // 交叉干活，时间点交叉，地块交叉，需要多次调试
        String did = "EC71BT2404140062";
        String yyyyMMdd = "20240426";
        String startTime = yyyyMMdd + "000000";
        String endTime = yyyyMMdd + "235959";
        double jobWidth = 3.5;
        测试拆分数据(did, startTime, endTime, jobWidth);
    }

    @Test
    void 测试1秒间隔015() {
        // 虽然路很多，但是都不符合地块作业的形态，所以也算不出来亩数
        String did = "NJ4GNBZAX0000160";
        String startTime = "20250603144925";
        String endTime = "20250603155017";
        double jobWidth = 2.4;
        测试拆分数据(did, startTime, endTime, jobWidth);
    }

    @Test
    void 测试1秒间隔016() {
        // 有一个地块由于他是作业关闭状态，所以没有识别出来，因为那个非作业状态上报时间间隔是5秒，一天中有两种上报时间间隔，只能计算上报时间间隔最短的那个
        String did = "NJ4GNBZAX0000273";
        String yyyyMMdd = "20250531";
        String startTime = yyyyMMdd + "000000";
        String endTime = yyyyMMdd + "235959";
        double jobWidth = 2.3;
        //测试拆分数据(did, startTime, endTime, jobWidth);
        测试拆分数据(did, startTime, endTime, jobWidth, false, new SplitRoadParams()
                .setDbScanEpsilon(18.0).setRoadWidth(3.0).setMinReturnMu(0.1).setCheckWorkingStatus(false));
    }

    @Test
    void 测试1秒间隔017() {
        // 识别出了3块地，有两块地中间因为路还符合了地块作业的形态，所以连上了
        String did = "EC71BT2404140062";
        String yyyyMMdd = "20240427";
        String startTime = yyyyMMdd + "000000";
        String endTime = yyyyMMdd + "235959";
        double jobWidth = 2.1;
        测试拆分数据(did, startTime, endTime, jobWidth);
    }

    @Test
    void 测试1秒间隔018() {
        // 无作业轨迹点
        String did = "EC71BT2404140062";
        String yyyyMMdd = "20240501";
        String startTime = yyyyMMdd + "000000";
        String endTime = yyyyMMdd + "235959";
        double jobWidth = 2.1;
        测试拆分数据(did, startTime, endTime, jobWidth);
    }

    @Test
    void 测试1秒间隔019() {
        // 识别出2块地
        String did = "EC71BT2404140062";
        String yyyyMMdd = "20240428";
        String startTime = yyyyMMdd + "000000";
        String endTime = yyyyMMdd + "235959";
        double jobWidth = 3.9;
        测试拆分数据(did, startTime, endTime, jobWidth);
    }

    @Test
    void 测试1秒间隔020() {
        // 识别出3块地
        String did = "EC71BT2404140062";
        String yyyyMMdd = "20240429";
        String startTime = yyyyMMdd + "000000";
        String endTime = yyyyMMdd + "235959";
        double jobWidth = 3.9;
        测试拆分数据(did, startTime, endTime, jobWidth);
    }

    @Test
    void 测试1秒间隔021() {
        // 识别出了多个地块，有一个地块由于非常细的一条，所以抛弃了，有一个地块由于点位密集度不够，断开了。
        String did = "JFT3352503S00207";
        String yyyyMMdd = "20251024";
        String startTime = yyyyMMdd + "000000";
        String endTime = yyyyMMdd + "235959";
        double jobWidth = 4;
        测试拆分数据(did, startTime, endTime, jobWidth, false, new SplitRoadParams().setDbScanEpsilon(12.0));
    }

    @Test
    void 测试1秒间隔022() {
        // 很多个地块，全都识别出来了，但是有一个地块，边缘少了一点点，多识别出来一块路
        String did = "EC73BD2504110767";
        String startTime = "20250529053827";
        String endTime = "20250529181257";
        double jobWidth = 2.6;
        测试拆分数据(did, startTime, endTime, jobWidth);
    }

    @Test
    void 测试1秒间隔023() {
        // 所有地块都识别出来了，但是有一块地连了一些路
        String did = "EC73BD2504110765";
        String startTime = "20250503055242";
        String endTime = "20250503181644";
        double jobWidth = 2.6;
        测试拆分数据(did, startTime, endTime, jobWidth);
    }

    @Test
    void 测试1秒间隔024() {
        // 识别出4块地，所有地块都能正确识别
        String did = "EC73BD2504110478";
        String startTime = "20250501071158";
        String endTime = "20250501183529";
        double jobWidth = 2.6;
        测试拆分数据(did, startTime, endTime, jobWidth);
    }

    @Test
    void 测试1秒间隔025() {
        // 所有地块都能识别出来
        String did = "EC73BD2504110477";
        String startTime = "20250501071207";
        String endTime = "20250501183550";
        double jobWidth = 2.6;
        测试拆分数据(did, startTime, endTime, jobWidth);
    }

    @Test
    void 测试1秒间隔026() {
        // 所有地块都能识别出来
        String did = "EC73BD2503190486";
        String startTime = "20250802055640";
        String endTime = "20250802164456";
        double jobWidth = 2.6;
        测试拆分数据(did, startTime, endTime, jobWidth);
    }

    @Test
    void 测试1秒间隔027() {
        // 所有地块都能识别出来
        String did = "EC73BD2504110122";
        String startTime = "20250522065331";
        String endTime = "20250522184733";
        double jobWidth = 2.6;
        测试拆分数据(did, startTime, endTime, jobWidth);
    }

    @Test
    void 测试1秒间隔028() {
        // 所有地块都能识别出来
        String did = "EC73BD2504110122";
        String startTime = "20250524115502";
        String endTime = "20250524220139";
        double jobWidth = 2.6;
        测试拆分数据(did, startTime, endTime, jobWidth);
    }

    @Test
    void 测试1秒间隔029() {
        // 所有地块都能识别出来
        String did = "EC73BD2504110774";
        String startTime = "20250526072025";
        String endTime = "20250526172848";
        double jobWidth = 2.5;
        测试拆分数据(did, startTime, endTime, jobWidth);
    }

    @Test
    void 测试1秒间隔030() {
        // 所有地块都能识别出来，但是带出来不少路上的亩数
        String did = "EC73BD2504110121";
        String startTime = "20250429063229";
        String endTime = "20250429175943";
        double jobWidth = 2.6;
        测试拆分数据(did, startTime, endTime, jobWidth, false, new SplitRoadParams().setRoadWidth(4.0));
    }

    @Test
    void 测试1秒间隔031() {
        // 所有地块都能识别出来
        String did = "EC73BD2503190486";
        String startTime = "20250421063426";
        String endTime = "20250421164529";
        double jobWidth = 2.6;
        测试拆分数据(did, startTime, endTime, jobWidth);
    }

    @Test
    void 测试1秒间隔032() {
        // 所有地块都能识别出来
        String did = "EC73BD2503190486";
        String startTime = "20250420065553";
        String endTime = "20250420175721";
        double jobWidth = 2.6;
        测试拆分数据(did, startTime, endTime, jobWidth);
    }

    @Test
    void 测试1秒间隔033() {
        // 所有地块都能识别出来
        String did = "EC71BD2501220049";
        String startTime = "20251031085228";
        String endTime = "20251031202509";
        double jobWidth = 2.8;
        测试拆分数据(did, startTime, endTime, jobWidth);
    }

    @Test
    void 测试1秒间隔034() {
        // 所有地块都能识别出来
        String did = "NJ4GNBZAX0000273";
        String startTime = "20250601000000";
        String endTime = "20250601235959";
        double jobWidth = 2.3;
        测试拆分数据(did, startTime, endTime, jobWidth);
    }

    @Test
    void 测试1秒间隔035() {
        // 所有地块都能识别出来
        String did = "EC73BD2504030112";
        String startTime = "20250414000000";
        String endTime = "20250414235959";
        double jobWidth = 2.6;
        测试拆分数据(did, startTime, endTime, jobWidth, false, new SplitRoadParams().setMinReturnMu(0.57));
    }

    @Test
    void 测试1秒间隔036() {
        // 所有地块都能识别出来，但多出来三条路也给算上亩数了
        String did = "EC73BD2504110478";
        String startTime = "20250527000000";
        String endTime = "20250527235959";
        double jobWidth = 2.6;
        测试拆分数据(did, startTime, endTime, jobWidth, false, new SplitRoadParams().setRoadWidth(4.0));
    }

    @Test
    void 测试1秒间隔037() {
        // 所有地块都能识别出来
        String did = "EC73BD2507220006";
        String startTime = "20251122064119";
        String endTime = "20251122172331";
        double jobWidth = 3;
        测试拆分数据(did, startTime, endTime, jobWidth, false);
    }

    @Test
    void 测试1秒间隔038() {
        // 所有地块都能识别出来
        String did = "EC73BD2507220006";
        String startTime = "20251123033002";
        String endTime = "20251123172703";
        double jobWidth = 3;
        测试拆分数据(did, startTime, endTime, jobWidth, false);
    }

    @Test
    void 测试1秒间隔039() {
        // 所有地块都能识别出来
        String did = "EC73BD2507220006";
        String startTime = "20251121065600";
        String endTime = "20251121175057";
        double jobWidth = 3;
        测试拆分数据(did, startTime, endTime, jobWidth, false);
    }

    @Test
    void 测试1秒间隔040() {
        // 所有地块都能识别出来
        String did = "EC73BD2507220006";
        String startTime = "20251120054614";
        String endTime = "20251120171321";
        double jobWidth = 3;
        测试拆分数据(did, startTime, endTime, jobWidth, false);
    }

    @Test
    void 测试1秒间隔041() {
        // 所有地块都能识别出来
        String did = "EC73BD2507220006";
        String startTime = "20251113094136";
        String endTime = "20251113164045";
        double jobWidth = 3;
        测试拆分数据(did, startTime, endTime, jobWidth, false);
    }

    @Test
    void 测试1秒间隔042() {
        // 所有地块都能识别出来
        String did = "EC73BD2507220006";
        String startTime = "20251110122141";
        String endTime = "20251110172259";
        double jobWidth = 3;
        测试拆分数据(did, startTime, endTime, jobWidth, false);
    }

    @Test
    void 测试1秒间隔043() {
        // 所有地块都能识别出来
        String did = "EC73BD2507220006";
        String startTime = "20251110075655";
        String endTime = "20251110113009";
        double jobWidth = 3;
        测试拆分数据(did, startTime, endTime, jobWidth, false);
    }

    @Test
    void 测试1秒间隔044() {
        // 所有地块都能识别出来
        String did = "EC73BD2507220006";
        String startTime = "20251122064119";
        String endTime = "20251122172331";
        double jobWidth = 3;
        测试拆分数据(did, startTime, endTime, jobWidth, false);
    }

    @Test
    void 测试1秒间隔045() {
        // 所有地块都能识别出来
        String did = "EC73BD2504110747";
        String startTime = "20251111000000";
        String endTime = "20251111235959";
        double jobWidth = 2.4;
        测试拆分数据(did, startTime, endTime, jobWidth);
    }


    @Test
    void 测试1秒间隔046() {
        // 所有地块都能识别出来
        String did = "EC73BD2504110136";
        String startTime = "20251103000000";
        String endTime = "20251103235959";
        double jobWidth = 2.4;
        测试拆分数据(did, startTime, endTime, jobWidth, false, new SplitRoadParams().setDbScanEpsilon(12.0));
    }

    @Test
    void 测试1秒间隔047() {
        // 所有地块都能识别出来
        String did = "EC73BD2509060268";
        String startTime = "20251111000000";
        String endTime = "20251111235959";
        double jobWidth = 2.4;
        测试拆分数据(did, startTime, endTime, jobWidth);
    }

    @Test
    void 测试1秒间隔048() {
        // 所有地块都能识别出来
        String did = "EC71BD2501220049";
        String startTime = "20251110072331";
        String endTime = "20251110190400";
        double jobWidth = 2.8;
        测试拆分数据(did, startTime, endTime, jobWidth, false);
    }

    @Test
    void 测试1秒间隔049() {
        // 这个是跨天的，测试一下
        String did = "EC71BD2501220049";
        String startTime = "20251114073955";
        String endTime = "20251115155306";
        double jobWidth = 2.8;
        测试拆分数据(did, startTime, endTime, jobWidth, false, new SplitRoadParams().setDbScanEpsilon(8.0));
    }


    @Test
    void 测试10秒间隔001() {
        // 一块地，非常多的路
        String did = "EM9101B8F5AZT0041";
        String yyyyMMdd = "20251027";
        String startTime = yyyyMMdd + "000000";
        String endTime = yyyyMMdd + "235959";
        double jobWidth = 3.5;
        测试拆分数据(did, startTime, endTime, jobWidth);
    }

    @Test
    void 测试10秒间隔002() {
        // 识别出4块地，有点漏地，给人家算少了
        String did = "EM9101B8F5AZT0041";
        String yyyyMMdd = "20251024";
        String startTime = yyyyMMdd + "000000";
        String endTime = yyyyMMdd + "235959";
        double jobWidth = 3.5;
        测试拆分数据(did, startTime, endTime, jobWidth);
    }

    @Test
    void 测试10秒间隔003() {
        // 形成不了作业地块，不是正常的作业轨迹
        String did = "EM9101B8F5AZT0041";
        String yyyyMMdd = "20250903";
        String startTime = yyyyMMdd + "000000";
        String endTime = yyyyMMdd + "235959";
        double jobWidth = 3.5;
        测试拆分数据(did, startTime, endTime, jobWidth);
    }

    @Test
    void 测试10秒间隔004() {
        // 识别出5块地，中间有个环形的地
        String did = "EM9101B8F5AZT0041";
        String yyyyMMdd = "20251025";
        String startTime = yyyyMMdd + "000000";
        String endTime = yyyyMMdd + "235959";
        double jobWidth = 3.5;
        测试拆分数据(did, startTime, endTime, jobWidth);
    }

    @Test
    void 测试10秒间隔005() {
        // 识别出3个地块
        String did = "EM9101B8F5AZT0041";
        String yyyyMMdd = "20251026";
        String startTime = yyyyMMdd + "000000";
        String endTime = yyyyMMdd + "235959";
        double jobWidth = 3.5;
        测试拆分数据(did, startTime, endTime, jobWidth);
    }

    @Test
    void 测试10秒间隔006() {
        // 识别出6个地块
        String did = "EM9101B8F5AZT0041";
        String yyyyMMdd = "20251023";
        String startTime = yyyyMMdd + "000000";
        String endTime = yyyyMMdd + "235959";
        double jobWidth = 2.5;
        测试拆分数据(did, startTime, endTime, jobWidth);
    }

    @Test
    void 测试10秒间隔007() {
        // 就跑了一圈就走了，不符合作业轨迹形态
        String did = "EM9101B8F4AZR0296";
        String yyyyMMdd = "20251101";
        String startTime = yyyyMMdd + "000000";
        String endTime = yyyyMMdd + "235959";
        double jobWidth = 3.5;
        测试拆分数据(did, startTime, endTime, jobWidth);
    }

    @Test
    void 测试15秒间隔001() {
        // 用户没有点击作业开始和结束，直接使用轨迹计算
        String did = "EC73BD2504030071";
        String startTime = "20250526000000";
        String endTime = "20250526235959";
        double jobWidth = 2.6;
        测试拆分数据(did, startTime, endTime, jobWidth, false, new SplitRoadParams().setCheckWorkingStatus(false));
    }

    @Test
    void 测试() {
        String did = "NJTEST0000000000";
        String startTime = "20260115000000";
        String endTime = "20260115235959";
        double jobWidth = 3.5;
        测试拆分数据(did, startTime, endTime, jobWidth);
    }

    @Test
    void 测试TdengineMapper() {
        TdengineMapper mapper = MyBatis.getMapper(TdengineMapper.class);
        List<DP> list = mapper.selectWorkPoints("EC73BD2509060248", LocalDateTimeUtil.parse("2025-10-23T00:00:00"), LocalDateTimeUtil.parse("2025-10-23T23:59:59"), true);
        log.info("{}", list.size());
    }

    @Test
    void 测试英轩农装设备1() {
        String did = "YXN26S2403T00066";
        String startTime = "20251223000000";
        String endTime = "20251223235959";
        double jobWidth = 3.5;
        测试拆分数据(did, startTime, endTime, jobWidth);
    }

    @Test
    void 测试大幅宽1() {
        String did = "EC73BD2504111716";
        String startTime = "20260115000000";
        String endTime = "20260115235959";
        double jobWidth = 20;
        测试拆分数据(did, startTime, endTime, jobWidth);
    }

    @Test
    void 测试大幅宽2() {
        String did = "EC73BD2504111716";
        String startTime = "20260115145741";
        String endTime = "20260115150038";
        double jobWidth = 6;
        List<DP> dps = selectWorkPoints(did, LocalDateTimeUtil.parse(startTime, "yyyyMMddHHmmss"), LocalDateTimeUtil.parse(endTime, "yyyyMMddHHmmss"));
        List<Wgs84Point> l = new ArrayList<>();
        for (DP dp : dps) {
            Wgs84Point wgs84Point = new Wgs84Point();
            wgs84Point.setGpsTime(dp.getP3014());
            wgs84Point.setLongitude(dp.getP2602());
            wgs84Point.setLatitude(dp.getP2603());
            if (dp.getP3020() != null) {
                wgs84Point.setJobStatus(2);//先设置非作业状态
                if (dp.getP3020() == 1) {// 终端ACC状态,0关闭，1开启
                    wgs84Point.setJobStatus(1);
                }
            }
            if (dp.getP4031() != null) {// 作业标识,1作业,0非作业,2暂停
                wgs84Point.setJobStatus(2);//先设置非作业状态
                if (dp.getP4031() == 1) {
                    wgs84Point.setJobStatus(1);//作业标识是1，认为是作业
                }
            }
            l.add(wgs84Point);
        }
        FarmPlot farmPlot = gisUtil.getFarmPlot(l, jobWidth);
        log.info("总点数 {}", farmPlot.getClusterPointCount());
        log.info("亩数 {}", farmPlot.getMu());
        log.info("作业时间 {} {}", farmPlot.getStartTime(), farmPlot.getEndTime());
        log.info("WKT: {}", farmPlot.getWkt());
        for (Wgs84Point wgs84Point : l) {
            System.out.println(StrUtil.format("{},{},{}", LocalDateTimeUtil.format(wgs84Point.getGpsTime(), "yyyyMMddHHmmss"), wgs84Point.getLongitude(), wgs84Point.getLatitude()));
        }
    }

    @Test
    void 测试getFarmPlot() {
        String did = "EC71BT2406060220";
        String startTime = "20251102130028";
        String endTime = "20251102153804";
        double jobWidth = 1.75;
        List<DP> dps = selectWorkPoints(did, LocalDateTimeUtil.parse(startTime, "yyyyMMddHHmmss"), LocalDateTimeUtil.parse(endTime, "yyyyMMddHHmmss"));
        List<Wgs84Point> l = new ArrayList<>();
        for (DP dp : dps) {
            Wgs84Point wgs84Point = new Wgs84Point();
            wgs84Point.setGpsTime(dp.getP3014());
            wgs84Point.setLongitude(dp.getP2602());
            wgs84Point.setLatitude(dp.getP2603());
            if (dp.getP3020() != null) {
                wgs84Point.setJobStatus(2);//先设置非作业状态
                if (dp.getP3020() == 1) {// 终端ACC状态,0关闭，1开启
                    wgs84Point.setJobStatus(1);
                }
            }
            if (dp.getP4031() != null) {// 作业标识,1作业,0非作业,2暂停
                wgs84Point.setJobStatus(2);//先设置非作业状态
                if (dp.getP4031() == 1) {
                    wgs84Point.setJobStatus(1);//作业标识是1，认为是作业
                }
            }
            l.add(wgs84Point);
        }
        FarmPlot farmPlot = gisUtil.getFarmPlot(l, jobWidth);
        log.info("总点数 {}", farmPlot.getClusterPointCount());
        log.info("亩数 {}", farmPlot.getMu());
        log.info("作业时间 {} {}", farmPlot.getStartTime(), farmPlot.getEndTime());
    }
}
