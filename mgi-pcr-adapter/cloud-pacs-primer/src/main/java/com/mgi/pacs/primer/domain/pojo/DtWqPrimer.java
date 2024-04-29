package com.mgi.pacs.primer.domain.pojo;

import com.alibaba.excel.annotation.ExcelProperty;
import com.baomidou.mybatisplus.annotation.IdType;
import com.baomidou.mybatisplus.annotation.TableField;
import com.baomidou.mybatisplus.annotation.TableId;
import com.baomidou.mybatisplus.annotation.TableName;
import io.swagger.annotations.ApiModel;
import io.swagger.annotations.ApiModelProperty;
import lombok.Data;
import lombok.EqualsAndHashCode;

import java.io.Serializable;
import java.util.Date;

/**
 * <p>
 * 引物表
 * </p>
 *
 * @author mabinbin
 * @since 2024-01-27
 */
@Data
@EqualsAndHashCode(callSuper = false)
@TableName("dt_wq_primer")
@ApiModel(value="DtWqPrimer对象", description="引物表")
public class DtWqPrimer implements Serializable {

    private static final long serialVersionUID = 1L;

    @ApiModelProperty(value = "主键id")
    @TableId(value = "id", type = IdType.AUTO)
    private Long id;

    @ApiModelProperty(value = "唯一对应一次实验")
    private String experiment_id;

    @ExcelProperty(value = "name_id")
    @ApiModelProperty(value = "化合物id")
    private String name_id;

    @ExcelProperty("chr_name1")
    @ApiModelProperty(value = "样本名称")
    private String chr_name1;

    @ExcelProperty("primer1")
    @ApiModelProperty(value = "碱基序列")
    private String primer1;

    @ExcelProperty("start1")
    @ApiModelProperty(value = "开始位置")
    private Integer start1;

    @ExcelProperty("end1")
    @ApiModelProperty(value = "结束位置")
    private Integer end1;

    @ExcelProperty("strand1")
    @ApiModelProperty(value = "序列方向")
    private String strand1;

    @ExcelProperty("chr_name2")
    @ApiModelProperty(value = "样本名称")
    private String chr_name2;

    @ExcelProperty("primer2")
    @ApiModelProperty(value = "碱基序列")
    private String primer2;

    @ExcelProperty("start2")
    @ApiModelProperty(value = "开始位置")
    private Integer start2;

    @ExcelProperty("end2")
    @ApiModelProperty(value = "结束位置")
    private Integer end2;

    @ExcelProperty("strand2")
    @ApiModelProperty(value = "序列方向")
    private String strand2;

    @ApiModelProperty(value = "创建时间")
    private Date create_time;

    @TableField(exist = false)
    @ExcelProperty("version")
    @ApiModelProperty(value = "版本号")
    private String version;
}
