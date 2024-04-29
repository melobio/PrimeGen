package com.mgi.pacs.common.core.feign;

import lombok.Data;
import org.springframework.beans.factory.annotation.Value;
import org.springframework.boot.context.properties.ConfigurationProperties;
import org.springframework.cloud.context.config.annotation.RefreshScope;
import org.springframework.context.annotation.Configuration;

import java.util.List;

/**
 * @author linzhunsheng
 * @date 2022/5/5
 */
@Data
@RefreshScope
@Configuration
@ConfigurationProperties("feign.inside")
public class FeignInsideAuthConfig {

    /**
     * feign请求前缀
     */
    public static final String FEIGN_INSIDE_URL_PREFIX = "/feign";

    /**
     * 内部调用 IP白名单
     */
    @Value("#{'${feign.inside.ips:}'.split(',')}")
    private List<String> ips;
}
