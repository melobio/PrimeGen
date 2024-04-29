package com.mgi.pacs.primer.domain;


import cn.hutool.core.util.ObjectUtil;
import com.alibaba.excel.annotation.ExcelProperty;
import com.mgi.pacs.primer.config.GenderConverter;
import com.mgi.pacs.primer.constant.CureKeysFromEnum;
import io.swagger.annotations.ApiModel;
import io.swagger.annotations.ApiModelProperty;
import lombok.Data;
import org.hibernate.validator.constraints.Range;

import javax.validation.constraints.NotBlank;
import javax.validation.constraints.NotNull;


@Data
@ApiModel("华大家人传过来的订单数据")
public class BgiFamilyOrderReq {

    @NotBlank(message = "订单号不能为空")
    @ApiModelProperty(value = "华大家人的订单号")
    private String orderId;

    @ExcelProperty(index = 0, value = "姓名")
    @NotBlank(message = "受检人姓名不能为空")
    @ApiModelProperty(value = "受检人姓名")
    private String name;

    @ExcelProperty(index = 1, value = "姓别", converter = GenderConverter.class)
    @Range(min = 0, max = 1, message = "性别只能取0-1")
    @ApiModelProperty(value = "受检人性别（0 女 1 男）")
    private Integer gender;

    @ExcelProperty(index = 2, value = "手机号")
    @NotBlank(message = "受检人电话不能为空")
    @ApiModelProperty(value = "受检人电话")
    private String phone;

    @NotNull(message = "年龄不能为空")
    @Range(min = 1, message = "年龄不能小于1岁")
    @ApiModelProperty(value = "年龄")
    private Integer age;

    @ExcelProperty(index = 3, value = "身份证")
    @ApiModelProperty(value = "身份证号")
    private String idCard;

    @NotBlank(message = "数据来源不能为空")
    @ApiModelProperty(value = "数据来源")
    private String from;

    @NotNull(message = "乳腺检查选项 0 不做 1 做")
    @Range(min = 0, max = 1, message = "乳腺检查选项 0 不做 1 做")
    @ApiModelProperty(value = "乳腺检查选项 0 不做 1 做")
    private Integer optionFF;

    @NotNull(message = "远程超声检查选项 0 不做 1 做")
    @Range(min = 0, max = 1, message = "远程超声检查选项 0 不做 1 做")
    @ApiModelProperty(value = "远程超声检查选项 0 不做 1 做")
    private Integer optionR3;

    @ApiModelProperty(value = "工号", example = "BGI123456")
    private String jobNum;

    //下面的传参目前不传
    @ApiModelProperty(value = "数据来源类型", hidden = true)
    private Integer fromType;

    @ApiModelProperty(value = "机构id")
    private Long institutionId;

    @ApiModelProperty(value = "社康id")
    private Long communityId;

    @ApiModelProperty(value = "套餐id")
    private Long setmealId;

    @ApiModelProperty(value = "R3 机构id")
    private Long institutionR3Id;

    @ApiModelProperty(value = "R3 社康id")
    private Long communityR3Id;

    @ApiModelProperty(value = "R3 套餐id")
    private Long setmealR3Id;

    //检查
    public void check() {

        //Assert.isTrue(PhoneUtil.isPhone(phone), () -> new CustomException("请输入正确的手机号"));
        this.fromType = CureKeysFromEnum.getTypeByName(from).getType();

        if (ObjectUtil.isEmpty(jobNum)) {
            this.jobNum = "";
        }
    }
}
