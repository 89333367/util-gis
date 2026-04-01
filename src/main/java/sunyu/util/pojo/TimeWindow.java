package sunyu.util.pojo;

import java.util.List;

/**
 * 时间窗口
 * <p>
 * 用于存储按时间间隔特征分割后的轨迹段
 */
public class TimeWindow {
    /**
     * 时间间隔（秒），表示该窗口内相邻点位的时间间隔特征
     */
    private long interval;
    /**
     * 该窗口内的WGS84点位集合
     */
    private List<Wgs84Point> points;

    public TimeWindow() {
    }

    /**
     * 构造方法
     *
     * @param interval 时间间隔（秒）
     * @param points   点位集合
     */
    public TimeWindow(long interval, List<Wgs84Point> points) {
        this.interval = interval;
        this.points = points;
    }

    public long getInterval() {
        return interval;
    }

    public void setInterval(long interval) {
        this.interval = interval;
    }

    public List<Wgs84Point> getPoints() {
        return points;
    }

    public void setPoints(List<Wgs84Point> points) {
        this.points = points;
    }
}
