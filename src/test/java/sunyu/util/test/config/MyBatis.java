package sunyu.util.test.config;

import cn.hutool.core.io.FileUtil;
import cn.hutool.core.lang.ClassScanner;
import cn.hutool.setting.dialect.Props;
import com.zaxxer.hikari.HikariConfig;
import com.zaxxer.hikari.HikariDataSource;
import org.apache.ibatis.builder.xml.XMLMapperBuilder;
import org.apache.ibatis.mapping.Environment;
import org.apache.ibatis.session.*;
import org.apache.ibatis.transaction.TransactionFactory;
import org.apache.ibatis.transaction.jdbc.JdbcTransactionFactory;

import javax.sql.DataSource;
import java.io.File;
import java.io.FileInputStream;
import java.lang.reflect.Proxy;
import java.net.URL;
import java.util.List;
import java.util.Set;
import java.util.concurrent.ConcurrentHashMap;

/**
 * 多数据源路由器
 *
 * @author SunYu
 */
public class MyBatis {
    private static final Props props = ConfigProperties.getProps();
    private static final ConcurrentHashMap<Class<?>, Object> MAPPER_CACHE = new ConcurrentHashMap<>();

    static {
        // 扫描 mapper 接口
        Set<Class<?>> mapperClasses = ClassScanner.scanPackageByAnnotation(props.getStr("mapper-package"), DS.class);
        for (Class<?> mapperClass : mapperClasses) {
            DS ds = mapperClass.getAnnotation(DS.class);
            if (ds != null) {
                registerMapper(mapperClass, ds.value());
            }
        }
    }


    private static void registerMapper(Class<?> mapperClass, String dsKey) {
        SqlSessionFactory factory = createFactory(dsKey, mapperClass);
        Object proxy = createProxy(mapperClass, factory);
        MAPPER_CACHE.put(mapperClass, proxy);
    }

    private static SqlSessionFactory createFactory(String dsKey, Class<?> mapperClass) {
        String prefix = dsKey + ".";

        // HikariCP 配置
        HikariConfig hikariConfig = new HikariConfig();
        hikariConfig.setDriverClassName(props.getStr(prefix + "driver"));
        hikariConfig.setJdbcUrl(props.getStr(prefix + "url"));
        hikariConfig.setUsername(props.getStr(prefix + "username"));
        hikariConfig.setPassword(props.getStr(prefix + "password"));
        hikariConfig.setMinimumIdle(0);
        hikariConfig.setMaximumPoolSize(10);

        // 创建数据源
        DataSource dataSource = new HikariDataSource(hikariConfig);

        // 创建事务工厂
        TransactionFactory transactionFactory = new JdbcTransactionFactory();

        // 创建环境
        Environment environment = new Environment(dsKey, transactionFactory, dataSource);

        Configuration configuration = new Configuration(environment);
        configuration.setMapUnderscoreToCamelCase(true);
        configuration.setDefaultExecutorType(ExecutorType.REUSE);
        configuration.setLogImpl(org.apache.ibatis.logging.slf4j.Slf4jImpl.class);

        // 加载 XML Mapper
        String xmlLocation = props.getStr(prefix + "mapperLocation");
        if (xmlLocation != null) {
            loadXmlMappers(configuration, xmlLocation);
        }

        // 只有在 XML 中没有定义该接口的映射时才添加接口
        // 避免重复注册导致的 BindingException
        if (!configuration.hasMapper(mapperClass)) {
            configuration.addMapper(mapperClass);
        }
        return new SqlSessionFactoryBuilder().build(configuration);
    }

    /**
     * 加载 XML Mapper
     */
    private static void loadXmlMappers(Configuration configuration, String locationPattern) {
        String baseDir = locationPattern.substring(0, locationPattern.lastIndexOf('/'));
        URL dirUrl = MyBatis.class.getClassLoader().getResource(baseDir);
        if (dirUrl == null) return;
        try {
            if ("file".equals(dirUrl.getProtocol())) {
                File dir = new File(dirUrl.toURI());
                List<File> xmlFiles = FileUtil.loopFiles(dir, file -> file.getName().endsWith(".xml"));

                for (File xmlFile : xmlFiles) {
                    String xmlPath = baseDir + "/" + xmlFile.getName();
                    try (FileInputStream is = new FileInputStream(xmlFile)) {
                        XMLMapperBuilder builder = new XMLMapperBuilder(
                                is, configuration, xmlPath, configuration.getSqlFragments()
                        );
                        builder.parse();
                    }
                }
            }
        } catch (Exception e) {
            throw new RuntimeException("加载 XML Mapper 失败: " + locationPattern, e);
        }
    }

    private static Object createProxy(Class<?> mapperClass, SqlSessionFactory factory) {
        return Proxy.newProxyInstance(
                mapperClass.getClassLoader(),
                new Class<?>[]{mapperClass},
                (proxy, method, args) -> {
                    // 特殊处理方法，避免空指针异常
                    if ("toString".equals(method.getName()) && method.getParameterCount() == 0) {
                        return mapperClass.getSimpleName() + "_Proxy";
                    }

                    try (SqlSession session = factory.openSession(true)) {
                        Object mapper = session.getMapper(mapperClass);
                        if (mapper == null) {
                            throw new RuntimeException("无法获取 Mapper 实例: " + mapperClass.getName());
                        }
                        return method.invoke(mapper, args);
                    } catch (Exception e) {
                        throw new RuntimeException("调用 Mapper 方法失败: " + method.getName(), e);
                    }
                }
        );
    }

    /**
     * 获取 Mapper 实例
     *
     * @param mapperClass Mapper 接口类
     * @param <T>         Mapper 接口类型
     * @return Mapper 实例
     */
    public static <T> T getMapper(Class<T> mapperClass) {
        return (T) MAPPER_CACHE.get(mapperClass);
    }
}