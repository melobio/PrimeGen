package com.mgi.pacs.primer.domain;

import io.swagger.annotations.ApiModel;
import io.swagger.annotations.ApiModelProperty;
import lombok.Data;

import javax.validation.constraints.Pattern;

@ApiModel("通用传参")
@Data
public class CommonBO {

    @Pattern(regexp = "^(?i)select\\s+[^;]+$", message = "SQL must start with 'select' and not contain ';'.")
    @ApiModelProperty("sql")
    private String sql;
}
