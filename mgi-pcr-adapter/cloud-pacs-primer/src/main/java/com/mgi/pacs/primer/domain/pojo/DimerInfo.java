package com.mgi.pacs.primer.domain.pojo;

import com.alibaba.excel.annotation.ExcelProperty;
import com.baomidou.mybatisplus.annotation.IdType;
import com.baomidou.mybatisplus.annotation.TableId;
import com.baomidou.mybatisplus.annotation.TableName;
import io.swagger.annotations.ApiModel;
import io.swagger.annotations.ApiModelProperty;
import lombok.Data;
import lombok.EqualsAndHashCode;

import java.io.Serializable;

/**
 * <p>
 * 二聚体信息
 * </p>
 *
 * @author mabinbin
 * @since 2024-01-27
 */
@Data
@EqualsAndHashCode(callSuper = false)
@TableName("dimer_info")
@ApiModel(value="DimerInfo对象", description="二聚体信息")
public class DimerInfo implements Serializable {

    private static final long serialVersionUID = 1L;

    @ApiModelProperty(value = "主键id")
    @TableId(value = "id", type = IdType.AUTO)
    private Long id;

    @ApiModelProperty(value = "唯一对应一次实验")
    private String experiment_id;

    @ApiModelProperty(value = "实验次序")
    private String experiment_order;

    @ExcelProperty("dimer")
    @ApiModelProperty(value = "二聚体序列")
    private String dimer;

    @ExcelProperty("r1_primer")
    @ApiModelProperty(value = "第一个引物名称")
    private String r1_primer;

    @ExcelProperty("r2_primer")
    @ApiModelProperty(value = "第二个引物名称")
    private String r2_primer;


}
