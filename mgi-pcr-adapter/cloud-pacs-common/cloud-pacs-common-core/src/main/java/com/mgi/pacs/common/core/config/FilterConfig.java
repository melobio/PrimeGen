package com.mgi.pacs.common.core.config;

import com.mgi.pacs.common.core.filter.IpWhiteFilter;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.boot.web.servlet.FilterRegistrationBean;
import org.springframework.context.annotation.Bean;
import org.springframework.context.annotation.Configuration;
import org.springframework.context.annotation.Lazy;

import javax.servlet.DispatcherType;

/**
 * 自定义过滤器配置类
 *
 * @author linzhunsheng
 * @date 2022/5/10
 */
@Configuration
public class FilterConfig {

    @Autowired
    private IpWhiteFilter ipWhiteFilter;

    @Bean
    @Lazy
    public FilterRegistrationBean<IpWhiteFilter> ipWhiteFilterRegistration() {
        FilterRegistrationBean<IpWhiteFilter> registration = new FilterRegistrationBean<>();
        // 添加过滤器
        registration.setFilter(ipWhiteFilter);
        // 设置过滤路径
        registration.addUrlPatterns("/feign/*", "/insider/*");
        registration.setName("ipWhiteFilter");
        // 设置优先级
        registration.setOrder(0);
        registration.setDispatcherTypes(DispatcherType.REQUEST);
        return registration;
    }
}
