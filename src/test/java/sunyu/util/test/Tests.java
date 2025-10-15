package sunyu.util.test;

import cn.hutool.core.convert.Convert;
import cn.hutool.core.date.DateUtil;
import cn.hutool.core.io.FileUtil;
import cn.hutool.core.map.MapUtil;
import sunyu.util.pojo.TrackPoint;

import java.util.ArrayList;
import java.util.List;
import java.util.TreeMap;

public class Tests {
    public static void main(String[] args) {
        List<TrackPoint> trajectory = getDatas(); // 包含经纬度、速度、方向、时间

    }

    private static List<TrackPoint> getDatas() {
        List<TrackPoint> l = new ArrayList<>();
        ProtocolSdk sdk = new ProtocolSdk("http://192.168.11.8/config.xml");
        for (String line : FileUtil.readUtf8Lines("d:/tmp/拆分算法测试/NJ4GNBZAX0000232/NJ4GNBZAX0000232_20250404.txt")) {
            TreeMap<String, String> m = sdk.parseProtocolString(line);
            if (MapUtil.isNotEmpty(m) && m.containsKey("2601") && m.get("2601").equals("0")) {
                TrackPoint t = new TrackPoint(DateUtil.parse(m.get("3014")).toJdkDate().getTime(),
                        Convert.toDouble(m.get("2602")), Convert.toDouble(m.get("2603")), Convert.toDouble(m.get("2204")), Convert.toInt(m.get("3012")));
                l.add(t);
            }
        }
        return l;
    }
}
