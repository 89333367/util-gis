package sunyu.util.test.config;

import cn.hutool.db.Entity;
import cn.hutool.db.handler.RsHandler;

import java.nio.charset.StandardCharsets;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.List;

public class TDengineHandler implements RsHandler<List<Entity>> {
    @Override
    public List<Entity> handle(ResultSet rs) throws SQLException {
        List<Entity> list = new ArrayList<>();
        while (rs.next()) {
            Entity entity = new Entity();
            int columnCount = rs.getMetaData().getColumnCount();

            for (int i = 1; i <= columnCount; i++) {
                String columnName = rs.getMetaData().getColumnLabel(i);
                Object value = rs.getObject(i);

                // 关键修复：TDengine VARCHAR字段被返回为byte[]的处理
                if (value instanceof byte[]) {
                    String columnTypeName = rs.getMetaData().getColumnTypeName(i);
                    // TDengine的字符串类型有VARCHAR, NCHAR, BINARY
                    if ("VARCHAR".equalsIgnoreCase(columnTypeName) ||
                            "NCHAR".equalsIgnoreCase(columnTypeName) ||
                            "BINARY".equalsIgnoreCase(columnTypeName)) {
                        value = new String((byte[]) value, StandardCharsets.UTF_8);
                    }
                }

                entity.set(columnName, value);
            }
            list.add(entity);
        }
        return list;
    }
}
