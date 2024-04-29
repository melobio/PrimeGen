package com.mgi.pacs.common.core.response;

/**
 * 响应枚举类
 *
 * @author linzhunsheng
 * @date 2022/03/16
 */
public enum ResponseEnum {

    /**
     * ok
     */
    OK(200, "ok"),

    /**
     * 一些需要登录的接口，而实际上因为前端无法知道token是否已过期，导致token已失效时，
     * 应该返回一个状态码，告诉前端token已经失效了，及时清理
     */
    CLEAN_TOKEN(10001, "clean token"),

    /**
     * 方法参数没有校验，内容由输入内容决定
     */
    METHOD_ARGUMENT_NOT_VALID(10002, ""),

    /**
     * 无法读取获取请求参数
     */
    HTTP_MESSAGE_NOT_READABLE(10003, "请求参数格式有误"),

    /**
     * 未授权
     */
    UNAUTHORIZED(10004, "Unauthorized"),

    /**
     * 服务器出了点小差
     */
    EXCEPTION(500, "服务器出了点小差"),

    /**
     * 用于直接显示提示用户的错误，内容由输入内容决定
     */
    SHOW_FAIL(10006, ""),

    /**
     * 数据异常
     */
    DATA_ERROR(10007, "数据异常，请刷新后重新操作"),

    /**
     * 数据不完整
     */
    DATA_INCOMPLETE(10008, "数据不完整"),

    /**
     * 没有查询权限
     */
    REFUND_NOT_PERMISSION(10009, "refund not permission"),

    /**
     * 云数据不存在
     */
    CLOUD_DATA_NOT_FOUND(10010, "cloud_data_not_found"),

    /**
     * 云数据状态更新失败
     */
    CLOUD_DATA_UPDATE_FAIL(10011, "cloud_data_update_fail");

    private final int code;

    private final String msg;

    public int value() {
        return code;
    }

    public String getMsg() {
        return msg;
    }

    ResponseEnum(int code, String msg) {
        this.code = code;
        this.msg = msg;
    }

    @Override
    public String toString() {
        return "ResponseEnum{" + "code='" + code + '\'' + ", msg='" + msg + '\'' + "} " + super.toString();
    }

}
