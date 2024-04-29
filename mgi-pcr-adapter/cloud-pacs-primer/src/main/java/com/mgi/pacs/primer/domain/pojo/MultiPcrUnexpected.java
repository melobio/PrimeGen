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
 * 多重PCR（Polymerase Chain Reaction，聚合酶链反应）实验中出现的意外或非预期结果
 * </p>
 *
 * @author mabinbin
 * @since 2024-01-27
 */
@Data
@EqualsAndHashCode(callSuper = false)
@TableName("multi_pcr_unexpected")
@ApiModel(value="MultiPcrUnexpected对象", description="多重PCR（Polymerase Chain Reaction，聚合酶链反应）实验中出现的意外或非预期结果")
public class MultiPcrUnexpected implements Serializable {

    private static final long serialVersionUID = 1L;

    @ApiModelProperty(value = "主键id")
    @TableId(value = "id", type = IdType.AUTO)
    private Long id;

    @ApiModelProperty(value = "唯一对应一次实验")
    private String experiment_id;

    @ApiModelProperty(value = "实验次序")
    private String experiment_order;

    @ExcelProperty("fprimer:rprimer")
    @ApiModelProperty(value = "聚合酶链反应（PCR）的两种特异性引物")
    private String fprimer_rprimer;

    @ExcelProperty("count")
    @ApiModelProperty(value = "引物数量")
    private Integer count;


}
