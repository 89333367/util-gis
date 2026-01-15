package sunyu.util.test.mapper.farm;

import org.apache.ibatis.annotations.Param;
import sunyu.util.test.config.DS;
import sunyu.util.test.entity.FarmWork;

import java.time.LocalDateTime;

@DS("farm")
public interface FarmMapper {

    void updateFarmWork(@Param("farmWork") FarmWork farmWork);

}
