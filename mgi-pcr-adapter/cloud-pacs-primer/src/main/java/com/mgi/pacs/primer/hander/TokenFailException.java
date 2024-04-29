package com.mgi.pacs.primer.hander;

/**
 * 自定义token异常类
 *
 * @author linzhunsheng
 * @date 2021/12/23
 */
public class TokenFailException extends RuntimeException {
    private static final long serialVersionUID = -8805986516272514755L;

    private Integer code;

    private String message;

    public TokenFailException(String message) {
        this.message = message;
    }

    public TokenFailException(String message, Integer code) {
        this.message = message;
        this.code = code;
    }

    public TokenFailException(String message, Throwable e) {
        super(message, e);
        this.message = message;
    }

    @Override
    public String getMessage() {
        return message;
    }

    public Integer getCode() {
        return code;
    }
}
