package sunyu.util.test;

import cn.hutool.core.codec.Base64;
import cn.hutool.core.date.DateUtil;
import cn.hutool.core.exceptions.UtilException;
import cn.hutool.core.map.MapUtil;
import cn.hutool.core.util.*;
import cn.hutool.cron.CronUtil;
import cn.hutool.http.HttpRequest;
import cn.hutool.log.Log;
import cn.hutool.log.LogFactory;
import org.w3c.dom.DOMException;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.NodeList;

import java.io.InputStream;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

/**
 * 内部协议工具
 */
public class ProtocolSdk {
    //key：内部协议数字key，value：协议配置
    public Map<String, Map<String, String>> idProtocolMap = new HashMap<>();
    //key：内部协议英文key，value：协议配置
    public Map<String, Map<String, String>> enProtocolMap = new HashMap<>();
    private final Log log = LogFactory.get();
    private String url = "http://192.168.11.8/config.xml";
    private InputStream inputStream;
    private final String DID_REG = "[a-zA-Z0-9]+[ ]?[a-zA-Z0-9]+";

    /**
     * 初始化，默认使用开发环境配置；http://192.168.11.8/config.xml
     */
    public ProtocolSdk() {
        log.info("初始化内部协议工具");
        init();
    }

    /**
     * 初始化
     *
     * @param inputStream
     */
    public ProtocolSdk(InputStream inputStream) {
        this.inputStream = inputStream;
        log.info("初始化内部协议工具");
        init();
    }

    /**
     * 初始化，使用自定义配置
     *
     * @param url 配置文件url
     */
    public ProtocolSdk(String url) {
        this.url = url;
        log.info("初始化内部协议工具 {}", url);
        init();
    }

    /**
     * 初始化，使用自定义配置，并且提供配置更新时间
     *
     * @param url
     * @param cron
     */
    public ProtocolSdk(String url, String cron) {
        this.url = url;
        log.info("初始化内部协议工具 {}", url);
        init();
        if (!CronUtil.getScheduler().isStarted()) {
            CronUtil.setMatchSecond(true);
            CronUtil.start();
        }
        CronUtil.schedule(cron, (Runnable) () -> {
            log.info("更新配置信息");
            init();
        });
    }

    /**
     * 初始化变量
     */
    private void init() {
        Map<String, Map<String, String>> tempIdProtocolMap = new HashMap<>();
        Map<String, Map<String, String>> tempEnProtocolMap = new HashMap<>();
        try {
            Document xml;
            if (inputStream != null) {
                xml = XmlUtil.readXML(inputStream);
            } else {
                xml = XmlUtil.readXML(HttpRequest.get(url).execute().bodyStream());
            }
            NodeList propList = XmlUtil.getNodeListByXPath("/Element/prop", xml);
            List<Element> props = XmlUtil.transElements(propList);
            for (Element prop : props) {
                Map<String, String> config = new HashMap<>();
                String id = prop.getAttribute("id");
                tempIdProtocolMap.put(id, config);
                List<Element> child = XmlUtil.transElements(prop.getChildNodes());
                for (Element element : child) {
                    config.put(element.getTagName(), element.getTextContent());
                }
            }

            idProtocolMap = tempIdProtocolMap;

            idProtocolMap.forEach((id, config) -> {
                String en = config.get("en");
                if (StrUtil.isNotBlank(en)) {
                    config.put("id", id);
                    tempEnProtocolMap.put(en, config);
                }
            });

            enProtocolMap = tempEnProtocolMap;
            log.info("配置信息加载完毕");
        } catch (UtilException e) {
            log.error("更新配置信息失败！");
            log.error(e);
        } catch (DOMException e) {
            log.error("更新配置信息失败！");
            log.error(e);
        }
    }

    /**
     * 解析值，一般用于前台展示，这里会根据convertResult把数字的值，转换成中文值，根据formatDateTime值，转换成格式化的值
     *
     * @param protocolParams 内部协议参数
     *
     * @return
     */
    public TreeMap<String, String> parseValue(Map<String, String> protocolParams) {
        TreeMap<String, String> m = new TreeMap<>();
        try {
            protocolParams.forEach((k, v) -> {
                Map<String, String> config = idProtocolMap.get(k);
                if (MapUtil.isEmpty(config)) {
                    config = idProtocolMap.get(protocolParams.get("params3") + "_" + k);
                }
                if (MapUtil.isEmpty(config)) {
                    config = enProtocolMap.get(k);
                }
                if (MapUtil.isNotEmpty(config)) {
                    String newValue = null;
                    if (config.getOrDefault("convertResult", "false").equals("true")) {//说明有需要转换的值
                        String ref = config.get("ref");//<ref>3003</ref>
                        if (StrUtil.isNotBlank(ref)) {//说明有引用
                            newValue = config.get("value_" + protocolParams.get(ref) + "_" + v);
                        } else {
                            String formatDateTime = config.get("formatDateTime");
                            if (StrUtil.isNotBlank(formatDateTime) && !v.contains(" // ")) {
                                try {
                                    newValue = DateUtil.parse(v).toString(formatDateTime);
                                } catch (Exception e) {
                                    log.error(e);
                                }
                            } else {
                                newValue = config.get("value_" + v);
                            }
                        }
                    }
                    if (StrUtil.isNotBlank(newValue)) {
                        v = newValue;
                    }
                }
                m.put(k, v);
            });
        } catch (Exception e) {
            log.error("解析失败，源：{}", protocolParams);
            log.error(e);
        }
        return m;
    }

    /**
     * 解析内部协议字符串；这里除了base64与base64Hex是解析过的，其余的信息返回值都是网关传过来的值，没有经过中文解析，一般用于存储hbase与hdfs
     * 这里还返回了id对应的英文key信息
     *
     * @param protocolString 内部协议字符串
     *
     * @return
     */
    public TreeMap<String, String> parseProtocolStringIncludeEnKey(String protocolString) {
        TreeMap<String, String> params = parseProtocolString(protocolString);
        if (MapUtil.isNotEmpty(params)) {
            TreeMap<String, String> newParams = new TreeMap<>();
            try {
                newParams.putAll(params);
                params.forEach((id, value) -> {
                    Map<String, String> config = idProtocolMap.get(id);
                    if (MapUtil.isEmpty(config)) {
                        config = idProtocolMap.get(params.get("params3") + "_" + id);
                    }
                    if (MapUtil.isNotEmpty(config)) {
                        if (StrUtil.isNotBlank(config.get("en"))) {
                            newParams.put(config.get("en"), value);
                        }
                    }
                });
            } catch (Exception e) {
                log.error("解析内部协议失败，源：{}", protocolString);
                log.error(e);
            }
            return newParams;
        }
        return params;
    }

    /**
     * 解析内部协议字符串；这里除了base64与base64Hex是解析过的，其余的信息返回值都是网关传过来的值，没有经过中文解析，一般用于存储hbase与hdfs
     *
     * @param protocolString 内部协议字符串
     *
     * @return
     */
    public TreeMap<String, String> parseProtocolString(String protocolString) {
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
                            value = convertBase64Value(params, id, value);//转换base64与base64Hex
                            params.put(id, value);// 使用内部协议的数字key
                        }
                    }
                }
                return params;
            }
        }
        return null;
    }

    private String convertBase64Value(Map<String, String> params, String id, String value) {
        Map<String, String> config = idProtocolMap.get(id);
        if (MapUtil.isEmpty(config)) {
            config = idProtocolMap.get(params.get("params3") + "_" + id);
        }
        if (MapUtil.isNotEmpty(config)) {
            if (config.getOrDefault("base64Hex", "false").equals("true")) {
                StringBuilder newValue = new StringBuilder();
                String[] vv = value.split("\\|");
                for (int i = 0; i < vv.length; i++) {
                    if (i > 0) {
                        newValue.append("|");
                    }
                    newValue.append(HexUtil.encodeHexStr(Base64.decode(vv[i])));
                }
                value = newValue.toString();
            } else if (config.getOrDefault("base64", "false").equals("true")) {
                StringBuilder newValue = new StringBuilder();
                String[] vv = value.split("\\|");
                for (int i = 0; i < vv.length; i++) {
                    if (i > 0) {
                        newValue.append("|");
                    }
                    newValue.append(Base64.decodeStr(vv[i], CharsetUtil.UTF_8));
                }
                value = newValue.toString();
            }
        }
        return value;
    }

    /**
     * 通过内部协议key获得中文描述
     *
     * @param key 可以是英文key或者数字key
     *
     * @return 如果没有找到中文描述，返回本身key
     */
    public String getCn(String key) {
        Map<String, String> protocolConfig;
        if (NumberUtil.isNumber(key)) {
            protocolConfig = idProtocolMap.get(key);
        } else {
            protocolConfig = enProtocolMap.get(key);
        }
        if (MapUtil.isNotEmpty(protocolConfig)) {
            String unit = protocolConfig.get("unit");
            if (StrUtil.isNotBlank(unit)) {
                return StrUtil.format("{} {}", protocolConfig.get("cn"), unit);
            } else {
                String cn = protocolConfig.get("cn");
                if (StrUtil.isNotBlank(cn)) {
                    return protocolConfig.get("cn");
                }
            }
        }
        return key;
    }

    /**
     * 通过英文key获得数字key
     *
     * @param en 英文key
     *
     * @return 数字key
     */
    public String getNum(String en) {
        Map<String, String> protocolConfig = enProtocolMap.get(en);
        if (MapUtil.isNotEmpty(protocolConfig)) {
            return protocolConfig.get("id");
        }
        return null;
    }

    /**
     * 通过数字key获得英文key
     *
     * @param num 数字key
     *
     * @return 英文key
     */
    public String getEn(String num) {
        Map<String, String> protocolConfig = idProtocolMap.get(num);
        if (MapUtil.isNotEmpty(protocolConfig)) {
            return protocolConfig.get("en");
        }
        return null;
    }

    /**
     * 通过数字key获得单位
     *
     * @param num 数字key
     *
     * @return 单位
     */
    public String getUnit(String num) {
        Map<String, String> protocolConfig = idProtocolMap.get(num);
        if (MapUtil.isNotEmpty(protocolConfig)) {
            return protocolConfig.get("unit");
        }
        return null;
    }

    /**
     * 将map转换成内部协议字符串
     *
     * @param datas
     *
     * @return
     */
    public String convertToProtocolString(Map<String, String> datas) {
        StringBuilder sb = new StringBuilder();
        // SUBMIT$6844369939890319360$WHX21137BK200409$REALTIME$TIME:20211003150957,gw:farm,3004:0
        sb.append(datas.get("params0"));
        sb.append("$");
        sb.append(datas.get("params1"));
        sb.append("$");
        sb.append(datas.get("params2"));
        sb.append("$");
        sb.append(datas.get("params3"));
        sb.append("$");
        sb.append("TIME:" + datas.get("TIME"));
        if (datas.containsKey("gw")) {
            sb.append(",gw:" + datas.get("gw"));
        }
        for (Map.Entry<String, String> e : datas.entrySet()) {
            if (NumberUtil.isNumber(e.getKey())) {
                sb.append("," + e.getKey() + ":" + e.getValue());
            }
        }
        return sb.toString();
    }

}
