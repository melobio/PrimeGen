package com.mgi.pacs.primer.config.swager;

import org.springframework.boot.autoconfigure.EnableAutoConfiguration;
import org.springframework.boot.autoconfigure.condition.ConditionalOnMissingBean;
import org.springframework.boot.autoconfigure.condition.ConditionalOnProperty;
import org.springframework.context.annotation.Bean;
import org.springframework.context.annotation.Configuration;
import springfox.documentation.builders.ApiInfoBuilder;
import springfox.documentation.builders.PathSelectors;
import springfox.documentation.builders.RequestHandlerSelectors;
import springfox.documentation.service.*;
import springfox.documentation.spi.DocumentationType;
import springfox.documentation.spi.service.contexts.SecurityContext;
import springfox.documentation.spring.web.plugins.ApiSelectorBuilder;
import springfox.documentation.spring.web.plugins.Docket;
import springfox.documentation.swagger2.annotations.EnableSwagger2;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.function.Predicate;

/**
 * @Description swagger配置类，需要配置自动配置
 **/

@Configuration
@EnableSwagger2 //开启swagger
@EnableAutoConfiguration
@ConditionalOnProperty(name = "swagger.enabled", matchIfMissing = true) //就算配置不存在也默认生效
public class SwaggerAutoConfiguration {

    /**
     * 默认的排除路径，排除Spring Boot默认的错误处理路径和端点
     */
    private static final List<String> DEFAULT_EXCLUDE_PATH = Arrays.asList("/error", "/actuator/**");

    private static final String BASE_PATH = "/**";


    //注入配置类
    @Bean
    @ConditionalOnMissingBean
    public SwaggerProperties swaggerProperties()
    {
        return new SwaggerProperties();
    }

    //编写api文档
    @Bean
    public Docket api(SwaggerProperties swaggerProperties) {

        // base-path处理
        if (swaggerProperties.getBasePath().isEmpty())
        {
            swaggerProperties.getBasePath().add(BASE_PATH);
        }

        // noinspection unchecked
        List<Predicate<String>> basePath = new ArrayList<Predicate<String>>();
        swaggerProperties.getBasePath().forEach(path -> basePath.add(PathSelectors.ant(path)));

        // exclude-path处理
        if (swaggerProperties.getExcludePath().isEmpty())
        {
            swaggerProperties.getExcludePath().addAll(DEFAULT_EXCLUDE_PATH);
        }

        List<Predicate<String>> excludePath = new ArrayList<>();
        swaggerProperties.getExcludePath().forEach(path -> excludePath.add(PathSelectors.ant(path)));

        //-------------------------------------------以上是扫描路径和不用扫描的路径---------------------------------------------------
        ApiSelectorBuilder builder = new Docket(DocumentationType.SWAGGER_2) //文档模板
                .host(swaggerProperties.getHost())  //自定义的host字符串
                .apiInfo(apiInfo(swaggerProperties)) //自定义apiInfo
                .select()
                .apis(RequestHandlerSelectors.basePackage(swaggerProperties.getBasePackage())); //扫描包的位置

        swaggerProperties.getBasePath().forEach(p -> builder.paths(PathSelectors.ant(p)));
        swaggerProperties.getExcludePath().forEach(p -> builder.paths(PathSelectors.ant(p).negate()));

        return builder
                .build()
                .securitySchemes(securitySchemes()) //构建请求头
                .securityContexts(securityContexts()) //安全上下文，全部接口都包含请求头
                .pathMapping("/");
    }

    /**
     * 自定义APIInfo
     * @param swaggerProperties 配置文件
     * @return apiInfo
     */
    private ApiInfo apiInfo(SwaggerProperties swaggerProperties)
    {
        return new ApiInfoBuilder()
                .title(swaggerProperties.getTitle())    //title
                .description(swaggerProperties.getDescription()) //描述
                .license(swaggerProperties.getLicense())    //licence
                .licenseUrl(swaggerProperties.getLicenseUrl())  //url
                .termsOfServiceUrl(swaggerProperties.getTermsOfServiceUrl()) //服务条款URL
                .contact(new Contact(swaggerProperties.getContact().getName(), swaggerProperties.getContact().getUrl(), swaggerProperties.getContact().getEmail()))
                .version(swaggerProperties.getVersion())    //版本
                .build();
    }

    /**
     * 安全模式，这里指定token通过Authorization头请求头传递
     */
    private List<SecurityScheme> securitySchemes()
    {
        List<SecurityScheme> apiKeyList = new ArrayList<SecurityScheme>();
        //apiKeyList.add(new ApiKey("Authorization", "Authorization", "header"));
        return apiKeyList;
    }

    /**
     * 安全上下文
     */
    private List<SecurityContext> securityContexts()
    {
        List<SecurityContext> securityContexts = new ArrayList<>();
        securityContexts.add(
                SecurityContext.builder()
                        .securityReferences(defaultAuth())
                        .operationSelector(o -> o.requestMappingPattern().matches("/.*"))
                        .build());
        return securityContexts;
    }

    /**
     * 默认的全局鉴权策略
     *
     * @return
     */
    private List<SecurityReference> defaultAuth()
    {
        AuthorizationScope authorizationScope = new AuthorizationScope("global", "accessEverything");
        AuthorizationScope[] authorizationScopes = new AuthorizationScope[1];
        authorizationScopes[0] = authorizationScope;

        List<SecurityReference> securityReferences = new ArrayList<>();
        //securityReferences.add(new SecurityReference("Authorization", authorizationScopes));
        return securityReferences;
    }

}
