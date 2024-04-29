package com.mgi.pacs.common.core.filter;

import cn.hutool.core.collection.CollectionUtil;
import cn.hutool.core.util.StrUtil;
import com.mgi.pacs.common.core.feign.FeignInsideAuthConfig;
import com.mgi.pacs.common.core.handler.HttpHandler;
import com.mgi.pacs.common.core.response.ResponseEnum;
import com.mgi.pacs.common.core.response.ServerResponseEntity;
import com.mgi.pacs.common.core.utils.IpUtils;
import lombok.extern.slf4j.Slf4j;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Component;

import javax.servlet.*;
import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;
import java.io.IOException;
import java.util.List;

/**
 * IP白名单过滤校验
 *
 * @author linzhunsheng
 * @date 2022/5/10
 */
@Slf4j
@Component
public class IpWhiteFilter implements Filter {

    @Autowired
    private FeignInsideAuthConfig feignInsideAuthConfig;

    @Autowired
    private HttpHandler httpHandler;

    @Override
    public void doFilter(ServletRequest request, ServletResponse response, FilterChain chain) throws IOException, ServletException {
        HttpServletRequest req = (HttpServletRequest) request;
        HttpServletResponse resp = (HttpServletResponse) response;

        // 拦截feign请求做IP校验
        if (!feignRequestCheck(req)) {
            httpHandler.printServerResponseToWeb(ServerResponseEntity.fail(ResponseEnum.UNAUTHORIZED));
            return;
        }

        chain.doFilter(req, resp);
    }

    private boolean feignRequestCheck(HttpServletRequest req) {
        // 不是feign请求，不用校验
        if (!req.getRequestURI().startsWith(FeignInsideAuthConfig.FEIGN_INSIDE_URL_PREFIX)) {
            return true;
        }

        // ip白名单
        List<String> ips = feignInsideAuthConfig.getIps();
        // 移除无用的空ip
        ips.removeIf(StrUtil::isBlank);
        // 有ip白名单，且ip不在白名单内，校验失败
        if (CollectionUtil.isNotEmpty(ips)
                && !ips.contains(IpUtils.getIpAddr())) {
            log.error("ip not in ip White list: {}, ip, {}", ips, IpUtils.getIpAddr());
            return false;
        }
        return true;
    }
}
