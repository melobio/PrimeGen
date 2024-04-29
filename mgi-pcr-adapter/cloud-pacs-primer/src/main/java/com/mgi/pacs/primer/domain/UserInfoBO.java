package com.mgi.pacs.primer.domain;

import lombok.Data;

import java.util.List;

@Data
public class UserInfoBO {

    private Integer code;
    private String msg;
    private DataBean data;

    @Data
    public static class DataBean {

        private Integer total;
        private List<RowsBean> rows;

        @Data
        public static class RowsBean {

            private String orderNo;
            private String name;
            private Integer gender;
            private String birthday;
            private Integer nation;
            private Integer country;
            private String telephone;
            private String idCard;
            private String beginTime;
            private String endTime;
            private Integer status;
            private String writeOffTime;
            private String creditCode;
            private String createTime;
        }
    }
}
