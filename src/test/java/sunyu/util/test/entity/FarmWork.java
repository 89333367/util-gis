package sunyu.util.test.entity;

import java.time.LocalDateTime;

public class FarmWork {
    private String did;
    private LocalDateTime jobEndTime;
    private Double effectiveJobArea;
    private String wktPoly;

    public String getDid() {
        return did;
    }

    public void setDid(String did) {
        this.did = did;
    }

    public LocalDateTime getJobEndTime() {
        return jobEndTime;
    }

    public void setJobEndTime(LocalDateTime jobEndTime) {
        this.jobEndTime = jobEndTime;
    }

    public Double getEffectiveJobArea() {
        return effectiveJobArea;
    }

    public void setEffectiveJobArea(Double effectiveJobArea) {
        this.effectiveJobArea = effectiveJobArea;
    }

    public String getWktPoly() {
        return wktPoly;
    }

    public void setWktPoly(String wktPoly) {
        this.wktPoly = wktPoly;
    }

}
