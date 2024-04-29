package com.mgi.pacs.common.core.page;

import io.swagger.annotations.ApiModelProperty;
import lombok.Data;

import java.util.List;

/**
 * 分页VO
 *
 * @author linzhunsheng
 * @date 2022/12/13
 */
@Data
public class PageVO<T> {

    @ApiModelProperty("总条目数")
    private Long total;

    @ApiModelProperty("结果集")
    private List<T> rows;
}
