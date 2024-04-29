package com.mgi.pacs.common.core.exception;

import com.mgi.pacs.common.core.response.ResponseEnum;

/**
 * 自定义异常类
 *
 * @author linzhunsheng
 * @date 2022/03/16
 */
public class CloudPacsException extends RuntimeException {

    private static final long serialVersionUID = 1L;

    private Object object;

    private ResponseEnum responseEnum;

    public CloudPacsException(String msg) {
        super(msg);
    }

    public CloudPacsException(String msg, Object object) {
        super(msg);
        this.object = object;
    }

    public CloudPacsException(String msg, Throwable cause) {
        super(msg, cause);
    }


    public CloudPacsException(ResponseEnum responseEnum) {
        super(responseEnum.getMsg());
        this.responseEnum = responseEnum;
    }

    public CloudPacsException(ResponseEnum responseEnum, Object object) {
        super(responseEnum.getMsg());
        this.responseEnum = responseEnum;
        this.object = object;
    }


    public Object getObject() {
        return object;
    }

    public ResponseEnum getResponseEnum() {
        return responseEnum;
    }

}
