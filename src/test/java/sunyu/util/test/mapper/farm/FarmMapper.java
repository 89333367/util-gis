package sunyu.util.test.mapper.farm;

import org.apache.ibatis.annotations.Param;
import sunyu.util.test.config.DS;
import sunyu.util.test.entity.FarmWork;
import sunyu.util.test.entity.FarmWorkSplitDay;

import java.util.List;

@DS("farm")
public interface FarmMapper {

    void updateFarmWork(@Param("farmWork") FarmWork farmWork);

    List<FarmWorkSplitDay> selectFarmWorkSplitDay(@Param("day") String day);
}
