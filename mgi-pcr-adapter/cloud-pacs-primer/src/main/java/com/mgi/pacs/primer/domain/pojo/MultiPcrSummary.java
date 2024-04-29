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
import java.math.BigDecimal;
import java.math.RoundingMode;

/**
 * <p>
 * 多重PCR统计数据
 * </p>
 *
 * @author mabinbin
 * @since 2024-01-27
 */
@Data
@EqualsAndHashCode(callSuper = false)
@TableName("multi_pcr_summary")
@ApiModel(value="MultipcrSummary对象", description="多重PCR统计数据")
public class MultiPcrSummary implements Serializable {

    private static final long serialVersionUID = 1L;

    @ApiModelProperty(value = "主键id")
    @TableId(value = "id", type = IdType.AUTO)
    private Long id;

    @ApiModelProperty(value = "唯一对应一次实验")
    private String experiment_id;

    @ExcelProperty("sample_name")
    @ApiModelProperty(value = "实验次数名称")
    private String sample_name;

    @ExcelProperty(value = "dimer_rate(%)")
    @ApiModelProperty(value = "二聚体的速率")
    private BigDecimal dimer_rate;

    @ExcelProperty(value = "uniformity(>0.1x)")
    @ApiModelProperty(value = "物质的均衡分布")
    private BigDecimal uniformity;

    @ExcelProperty(value = "map_target_rate(%)")
    @ApiModelProperty(value = "基因表达映射")
    private BigDecimal map_target_rate;

    @ExcelProperty(value = "cov_100x(%)")
    @ApiModelProperty(value = "数据质量评估")
    private BigDecimal cov_100x;
}
