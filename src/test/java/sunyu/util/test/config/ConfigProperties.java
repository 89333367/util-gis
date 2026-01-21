package sunyu.util.test.config;

import cn.hutool.core.util.StrUtil;
import cn.hutool.log.Log;
import cn.hutool.log.LogFactory;
import cn.hutool.setting.dialect.Props;

public class ConfigProperties {
    static Log log = LogFactory.get();
    static Props props;
    static String CONFIG_FILE_NAME = "application";

    static {
        log.info("读取 properties 文件");
        props = new Props(StrUtil.format("{}.properties", CONFIG_FILE_NAME));
        String env = props.getStr("env");
        if (StrUtil.isNotBlank(env)) {
            log.info("合并 {} 环境配置", env);
            props.putAll(new Props(StrUtil.format("{}-{}.properties", CONFIG_FILE_NAME, env)));
        }
        log.info("配置 {}", props);
    }

    public static Props getProps() {
        return props;
    }
}
