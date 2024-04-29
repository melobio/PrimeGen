package com.mgi.pacs.primer.hander;

import cn.hutool.core.exceptions.ValidateException;
import com.mgi.pacs.common.core.exception.CloudPacsException;
import com.mgi.pacs.common.core.response.ServerResponseEntity;
import com.mgi.pacs.primer.exception.BaseException;
import lombok.extern.slf4j.Slf4j;
import org.apache.commons.lang3.StringUtils;
import org.apache.shiro.authc.IncorrectCredentialsException;
import org.apache.shiro.authc.UnknownAccountException;
import org.apache.shiro.authz.AuthorizationException;
import org.apache.shiro.authz.UnauthorizedException;
import org.hibernate.validator.internal.engine.path.PathImpl;
import org.mybatis.spring.MyBatisSystemException;
import org.springframework.dao.DataIntegrityViolationException;
import org.springframework.http.converter.HttpMessageNotReadableException;
import org.springframework.web.HttpMediaTypeNotSupportedException;
import org.springframework.web.HttpRequestMethodNotSupportedException;
import org.springframework.web.bind.MethodArgumentNotValidException;
import org.springframework.web.bind.MissingServletRequestParameterException;
import org.springframework.web.bind.annotation.ExceptionHandler;
import org.springframework.web.bind.annotation.ResponseBody;
import org.springframework.web.bind.annotation.RestControllerAdvice;
import org.springframework.web.multipart.MultipartException;

import javax.servlet.http.HttpServletRequest;
import javax.validation.ConstraintViolation;
import javax.validation.ConstraintViolationException;
import javax.validation.UnexpectedTypeException;
import java.sql.SQLSyntaxErrorException;
import java.util.Comparator;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

/**
 * 会诊全局异常处理器
 *
 */
@Slf4j
@RestControllerAdvice
public class GlobalExceptionHandler {

    @ExceptionHandler(value = Exception.class)
    public ServerResponseEntity errorHandler(Exception ex) {
        log.error(ex.getMessage(), ex);
        if (ex instanceof BaseException) {
            return ServerResponseEntity .fail(ex.getMessage());
        }
        if (ex instanceof IllegalArgumentException) {
            return ServerResponseEntity.fail(ex.getMessage());
        }
        return ServerResponseEntity.fail(ex.getMessage());
    }

    @ExceptionHandler(value = SQLSyntaxErrorException.class)
    public ServerResponseEntity errorHandler(SQLSyntaxErrorException ex) {
        log.error(ex.getMessage(), ex);
        return ServerResponseEntity.fail("SQL出错" + ex.getMessage());
    }

    @ExceptionHandler(value = MyBatisSystemException.class)
    public ServerResponseEntity errorHandler(MyBatisSystemException ex) {
        log.error(ex.getMessage(), ex);
        return ServerResponseEntity.fail("持久层出错" + ex.getMessage());
    }

    /**
     * 上传文件异常
     */
    @ExceptionHandler(MultipartException.class)
    public Object validExceptionHandler(MultipartException e)
    {
        log.error(e.getMessage(), e);
        return ServerResponseEntity.fail(e.getMessage());
    }

    @ExceptionHandler(value = DataIntegrityViolationException.class)
    public ServerResponseEntity errorHandler(DataIntegrityViolationException ex) {
        log.error(ex.getMessage(), ex);
        return ServerResponseEntity.fail("持久层出错 字段数据太长" + ex.getMessage());
    }

    @ExceptionHandler(value = HttpMediaTypeNotSupportedException.class)
    public ServerResponseEntity errorHandler(HttpMediaTypeNotSupportedException ex) {
        log.error(ex.getMessage(), ex);
        return ServerResponseEntity.fail("请求数据类型不支持" + ex.getMessage());
    }

    @ResponseBody
    @ExceptionHandler(value = MissingServletRequestParameterException.class)
    public ServerResponseEntity errorHandler(MissingServletRequestParameterException ex) {
        log.error(ex.getMessage(), ex);
        return ServerResponseEntity.fail("参数缺失" + ex.getMessage());
    }

    @ResponseBody
    @ExceptionHandler(value = HttpMessageNotReadableException.class)
    public ServerResponseEntity errorHandler(HttpMessageNotReadableException ex) {
        log.error(ex.getMessage(), ex);
        return ServerResponseEntity.fail("请求数据格式错误" + ex.getMessage());
    }

    @ExceptionHandler(value = UnknownAccountException.class)
    public ServerResponseEntity errorHandler(UnknownAccountException ex) {
        log.error(ex.getMessage(), ex);
        return ServerResponseEntity.fail(ex.getMessage());
    }

    @ExceptionHandler(value = IncorrectCredentialsException.class)
    public ServerResponseEntity errorHandler(IncorrectCredentialsException ex) {
        log.error(ex.getMessage(), ex);
        return ServerResponseEntity.fail(ex.getMessage());
    }

    @ExceptionHandler(value = HttpRequestMethodNotSupportedException.class)
    public ServerResponseEntity errorHandler(HttpServletRequest request, HttpRequestMethodNotSupportedException ex) {
        log.error("request地址{}{}{}", request.getRequestURI(), ex.getMessage(), ex);
        return ServerResponseEntity.fail(ex.getMessage());
    }

    /**
     * 自定义验证异常
     */
    @ExceptionHandler(MethodArgumentNotValidException.class)
    public Object validExceptionHandler(MethodArgumentNotValidException e) {
        log.error(e.getMessage(), e);
        String message = e.getBindingResult().getFieldError().getDefaultMessage();
        return ServerResponseEntity.fail(message);
    }


    /**
     * 自定义验证异常
     */
    @ExceptionHandler(UnexpectedTypeException.class)
    public Object validExceptionHandler(UnexpectedTypeException e) {
        log.error(e.getMessage(), e);
        String message = e.getMessage();
        return ServerResponseEntity.fail(message);
    }

    @ExceptionHandler(value = ConstraintViolationException.class)
    public ServerResponseEntity handleMethodConstraintViolationException(ConstraintViolationException ex) {
        log.error(ex.getMessage());
        Set<ConstraintViolation<?>> constraintViolations = ex.getConstraintViolations();
        //排序 为了不让有多个报错的时候,每次返回的内容顺序不一致.
        List<ConstraintViolation> collect = constraintViolations.stream()
                .sorted(Comparator.comparing(ConstraintViolation::getMessage)).collect(Collectors.toList());
        for (ConstraintViolation<?> constraintViolation : constraintViolations) {
            PathImpl pathImpl = (PathImpl) constraintViolation.getPropertyPath();
            // 读取参数字段，constraintViolation.getMessage() 读取验证注解中的message值
            String paramName = pathImpl.getLeafNode().getName();
            String message = "错误信息为：".concat(constraintViolation.getMessage());
            log.error("{} -> {} -> {}", constraintViolation.getRootBeanClass().getName(), pathImpl.toString(), message);
            return ServerResponseEntity.fail(message);
        }
        return ServerResponseEntity.fail(ex.getMessage());
    }

    /**
     * 基础异常
     */
    @ExceptionHandler(BaseException.class)
    public ServerResponseEntity baseException(BaseException e) {
        return ServerResponseEntity.fail(e.getDefaultMessage());
    }

    /**
     * 参数校验异常
     */
    @ExceptionHandler(ValidateException.class)
    public ServerResponseEntity validateExceptionHandler(TokenFailException e) {
        return ServerResponseEntity.fail(e.getMessage());
    }
}
