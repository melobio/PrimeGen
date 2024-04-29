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

/**
 * <p>
 * 多重PCR(聚合酶链反应)实验中的效率
 * </p>
 *
 * @author mabinbin
 * @since 2024-01-27
 */
@Data
@EqualsAndHashCode(callSuper = false)
@TableName("multi_pcr_efficiency")
@ApiModel(value="MultiPcrEfficiency对象", description="多重PCR(聚合酶链反应)实验中的效率")
public class MultiPcrEfficiency implements Serializable {

    private static final long serialVersionUID = 1L;

    @ApiModelProperty(value = "主键id")
    @TableId(value = "id", type = IdType.AUTO)
    private Long id;

    @ApiModelProperty(value = "唯一对应一次实验")
    private String experiment_id;

    @ApiModelProperty(value = "实验次序")
    private String experiment_order;

    @ExcelProperty("amp_name")
    @ApiModelProperty(value = "序列名称")
    private String amp_name;

    @ExcelProperty("f_bias")
    @ApiModelProperty(value = "PCR反应中使用的正向引物的偏好性")
    private Integer f_bias;

    @ExcelProperty("r_bias")
    @ApiModelProperty(value = "反向引物在PCR反应中的偏好性")
    private Integer r_bias;

    @ExcelProperty("primer_efficiency(%)")
    @ApiModelProperty(value = "引物效率")
    private BigDecimal primer_efficiency;

}
