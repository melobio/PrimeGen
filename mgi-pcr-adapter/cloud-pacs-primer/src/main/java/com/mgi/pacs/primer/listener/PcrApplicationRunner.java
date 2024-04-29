package com.mgi.pacs.adapter.listenler;

import com.google.protobuf.ServiceException;
import lombok.SneakyThrows;
import lombok.extern.slf4j.Slf4j;
import org.apache.ibatis.io.Resources;
import org.apache.ibatis.jdbc.ScriptRunner;
import org.springframework.context.ApplicationContextInitializer;
import org.springframework.context.ConfigurableApplicationContext;

import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.SQLException;
import java.sql.Statement;

/*
 * @description: 创建数据库
 */

@Slf4j
public class PcrApplicationRunner implements ApplicationContextInitializer {

    @SneakyThrows
    @Override
    public void initialize(ConfigurableApplicationContext applicationContext) {

        Connection conn = getConnection();
        Statement stmt = null;
        try {
            conn.setAutoCommit(false);
            stmt = conn.createStatement();

            String createDatabaseSQL = "CREATE DATABASE IF NOT EXISTS `" + "mgi_primer" + "` DEFAULT CHARSET utf8mb4 COLLATE utf8mb4_general_ci;";
            stmt.execute(createDatabaseSQL);

            conn.setCatalog("mgi_primer");

            // 获取当前数据库名称
            log.info("当前数据库：{}", conn.getCatalog()); // 若未选择数据库，则 getCatalog 返回空

            //执行初始化脚本
            ScriptRunner runner = new ScriptRunner(conn);
            runner.setErrorLogWriter(null);
            runner.setLogWriter(null);
            runner.runScript(Resources.getResourceAsReader("init-sql-script/mgi_primer_init.sql"));
            conn.commit();
        } catch (Exception e) {
            e.printStackTrace();
            conn.rollback();
            throw new ServiceException("初始化用户数据脚本时出错");
        } finally {
            if (stmt != null) {
                stmt.close();
            }
            if (conn != null) {
                conn.close();
            }
        }
    }

    /**
     * 获取数据库连接
     */
    private Connection getConnection() throws SQLException {

        String SPRING_PCR_MYSQL_URL = System.getenv("SPRING_PCR_MYSQL_URL");
        String SPRING_PCR_MYSQL_USERNAME = System.getenv("SPRING_PCR_MYSQL_USERNAME");
        String SPRING_PCR_MYSQL_PASSWORD = System.getenv("SPRING_PCR_MYSQL_PASSWORD");

        return DriverManager.getConnection(SPRING_PCR_MYSQL_URL, SPRING_PCR_MYSQL_USERNAME, SPRING_PCR_MYSQL_PASSWORD);
    }
}
