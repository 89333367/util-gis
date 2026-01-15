package sunyu.util.test.mapper.tdengine;

import org.apache.ibatis.annotations.Param;
import sunyu.util.test.config.DS;
import sunyu.util.test.entity.DP;

import java.time.LocalDateTime;
import java.util.List;

@DS("tdengine")
public interface TdengineMapper {

    List<DP> selectWorkPoints(@Param("did") String did, @Param("startTime") LocalDateTime startTime, @Param("endTime") LocalDateTime endTime, @Param("checkWorkStatus") Boolean checkWorkStatus);
    List<DP> selectWorkProtocol(@Param("did") String did, @Param("startTime") LocalDateTime startTime, @Param("endTime") LocalDateTime endTime, @Param("checkWorkStatus") Boolean checkWorkStatus);

}
