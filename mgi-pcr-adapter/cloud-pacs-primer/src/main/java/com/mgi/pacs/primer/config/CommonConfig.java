package com.mgi.pacs.primer.config;

import org.mybatis.spring.annotation.MapperScan;
import org.springframework.context.annotation.Configuration;

@MapperScan(basePackages = "com.mgi.pacs.primer.mapper")
@Configuration
public class CommonConfig {

}
