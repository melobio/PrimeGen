package com.mgi.pacs.primer.domain.pojo;

import com.baomidou.mybatisplus.annotation.IdType;
import com.baomidou.mybatisplus.annotation.TableId;
import com.baomidou.mybatisplus.annotation.TableName;
import io.swagger.annotations.ApiModel;
import io.swagger.annotations.ApiModelProperty;
import lombok.Data;
import lombok.EqualsAndHashCode;

import java.io.Serializable;
import java.math.BigDecimal;

/**
 * <p>
 * 多重PCR的深度
 * </p>
 *
 * @author mabinbin
 * @since 2024-01-27
 */
@Data
@EqualsAndHashCode(callSuper = false)
@TableName("multi_pcr_depth")
@ApiModel(value="MultipcrDepth对象", description="多重PCR的深度")
public class MultiPcrDepth implements Serializable {

    private static final long serialVersionUID = 1L;

    @ApiModelProperty(value = "主键id")
    @TableId(value = "id", type = IdType.AUTO)
    private Long id;

    @ApiModelProperty(value = "唯一对应一次实验")
    private String experiment_id;

    @ApiModelProperty(value = "实验次序")
    private String experiment_order;

    @ApiModelProperty(value = "二聚物名称")
    private String name;

    @ApiModelProperty(value = "测序深度")
    private BigDecimal depth;


}
