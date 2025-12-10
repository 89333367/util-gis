package sunyu.util.pojo;

import java.util.ArrayList;
import java.util.List;

/**
 * 道路拆分结果
 */
public class SplitResult {
    /**
     * 轮廓的WKT表示（WGS84坐标系）
     */
    private String wkt;
    /**
     * 使用作业总宽幅（米）
     */
    private double workingWidth;

    /**
     * 整体轮廓面积（亩）
     */
    private double mu;

    /**
     * 拆分后的地块列表
     */
    private List<Part> parts = new ArrayList<>();

    public String getWkt() {
        return wkt;
    }

    public void setWkt(String wkt) {
        this.wkt = wkt;
    }

    public double getWorkingWidth() {
        return workingWidth;
    }

    public void setWorkingWidth(double workingWidth) {
        this.workingWidth = workingWidth;
    }

    public double getMu() {
        return mu;
    }

    public void setMu(double mu) {
        this.mu = mu;
    }

    public List<Part> getParts() {
        return parts;
    }

    public void setParts(List<Part> parts) {
        this.parts = parts;
    }
}
