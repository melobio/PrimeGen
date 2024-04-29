package com.mgi.pacs.common.core.handler;

import cn.hutool.core.util.StrUtil;
import com.mgi.pacs.common.core.exception.CloudPacsException;
import com.mgi.pacs.common.core.response.ResponseEnum;
import com.mgi.pacs.common.core.response.ServerResponseEntity;
import io.seata.core.context.RootContext;
import io.seata.core.exception.TransactionException;
import io.seata.tm.api.GlobalTransactionContext;
import lombok.extern.slf4j.Slf4j;
import org.hibernate.validator.internal.engine.path.PathImpl;
import org.springframework.http.HttpStatus;
import org.springframework.http.ResponseEntity;
import org.springframework.http.converter.HttpMessageNotReadableException;
import org.springframework.validation.BindException;
import org.springframework.validation.FieldError;
import org.springframework.web.bind.MethodArgumentNotValidException;
import org.springframework.web.bind.MissingServletRequestParameterException;
import org.springframework.web.bind.annotation.ExceptionHandler;
import org.springframework.web.bind.annotation.RestControllerAdvice;

import javax.validation.ConstraintViolation;
import javax.validation.ConstraintViolationException;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

/**
 * 自定义错误处理器，除了授权只要请求服务器成功，就返回200，错误根据错误码前端进行处理
 *
 * @author linzhunsheng
 * @date 2022/3/16
 */
@Slf4j
@RestControllerAdvice
public class DefaultExceptionHandlerConfig {

	@ExceptionHandler({ MethodArgumentNotValidException.class, BindException.class })
	public ResponseEntity<ServerResponseEntity<List<String>>> methodArgumentNotValidExceptionHandler(Exception e) {
		log.error("methodArgumentNotValidExceptionHandler", e);
		List<FieldError> fieldErrors = null;
		if (e instanceof MethodArgumentNotValidException) {
			fieldErrors = ((MethodArgumentNotValidException) e).getBindingResult().getFieldErrors();
		}
		if (e instanceof BindException) {
			fieldErrors = ((BindException) e).getBindingResult().getFieldErrors();
		}
		if (fieldErrors == null) {
			return ResponseEntity.status(HttpStatus.OK)
					.body(ServerResponseEntity.fail(ResponseEnum.METHOD_ARGUMENT_NOT_VALID));
		}

		List<String> defaultMessages = new ArrayList<>(fieldErrors.size());
		for (FieldError fieldError : fieldErrors) {
			defaultMessages.add(fieldError.getField() + ":" + fieldError.getDefaultMessage());
		}
		return ResponseEntity.status(HttpStatus.OK)
				.body(ServerResponseEntity.fail(ResponseEnum.METHOD_ARGUMENT_NOT_VALID, defaultMessages));
	}

	@ExceptionHandler(value = MissingServletRequestParameterException.class)
	public ResponseEntity<ServerResponseEntity<List<String>>> missingServletRequestParameterExceptionHandler(Exception e) {
		log.error("missingServletRequestParameterExceptionHandler", e);
		String parameterName = ((MissingServletRequestParameterException) e).getParameterName();
		return ResponseEntity.status(HttpStatus.OK)
				.body(ServerResponseEntity.fail("参数缺失：" + parameterName));
	}

	@ExceptionHandler(value = ConstraintViolationException.class)
	public ResponseEntity<ServerResponseEntity<Object>> methodConstraintViolationExceptionHandler(Exception ex) {
		Set<ConstraintViolation<?>> constraintViolations = ((ConstraintViolationException) ex).getConstraintViolations();
		List<ConstraintViolation<?>> collect = constraintViolations.stream()
				.sorted(Comparator.comparing(ConstraintViolation::getMessage)).collect(Collectors.toList());
		for (ConstraintViolation<?> constraintViolation : collect) {
			PathImpl pathImpl = (PathImpl) constraintViolation.getPropertyPath();
			String message = constraintViolation.getMessage();
			log.error("{} -> {} -> {}", constraintViolation.getRootBeanClass().getName(), pathImpl.toString(), message);
			return ResponseEntity.status(HttpStatus.OK).body(ServerResponseEntity.fail(message));
		}
		return ResponseEntity.status(HttpStatus.OK)
				.body(ServerResponseEntity.fail(ex.getMessage()));
	}

	@ExceptionHandler({ HttpMessageNotReadableException.class })
	public ResponseEntity<ServerResponseEntity<List<FieldError>>> methodArgumentNotValidExceptionHandler(
			HttpMessageNotReadableException e) {
		log.error("methodArgumentNotValidExceptionHandler", e);
		return ResponseEntity.status(HttpStatus.OK)
				.body(ServerResponseEntity.fail(ResponseEnum.HTTP_MESSAGE_NOT_READABLE));
	}

	/**
	 * pacs错误码
	 */
	@ExceptionHandler(CloudPacsException.class)
	public ResponseEntity<ServerResponseEntity<Object>> defaultExceptionHandler(CloudPacsException e) {
		log.error("DefaultExceptionHandler", e);

		ResponseEnum responseEnum = e.getResponseEnum();
		// 失败返回失败消息 + 状态码
		if (responseEnum != null) {
			return ResponseEntity.status(HttpStatus.OK).body(ServerResponseEntity.fail(responseEnum, e.getObject()));
		}
		// 失败返回消息 状态码固定为直接显示消息的状态码
		return ResponseEntity.status(HttpStatus.OK).body(ServerResponseEntity.showFailMsg(e.getMessage()));
	}

	@ExceptionHandler(Exception.class)
	public ResponseEntity<ServerResponseEntity<Object>> exceptionHandler(Exception e) throws TransactionException {
		log.error("exceptionHandler", e);
		log.info("RootContext.getXID(): " + RootContext.getXID());
		if (StrUtil.isNotBlank(RootContext.getXID())) {
			GlobalTransactionContext.reload(RootContext.getXID()).rollback();
		}
		return ResponseEntity.status(HttpStatus.OK).body(ServerResponseEntity.fail(ResponseEnum.EXCEPTION));
	}
}
