package com.mgi.pacs.common.datasource.config;

import org.mybatis.spring.annotation.MapperScan;
import org.springframework.context.annotation.Configuration;
import org.springframework.transaction.annotation.EnableTransactionManagement;

/**
 * mybatis-plus配置类
 *
 * @author linzhunsheng
 * @date 2022/01/10
 */
@Configuration
@MapperScan(basePackages = {"com.mgi.**.mapper"})
@EnableTransactionManagement
public class MybatisPlusConfig {

}
