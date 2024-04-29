package com.mgi.pacs.primer.config;

import cn.hutool.core.util.StrUtil;
import com.alibaba.excel.converters.Converter;
import com.alibaba.excel.metadata.GlobalConfiguration;
import com.alibaba.excel.metadata.data.ReadCellData;
import com.alibaba.excel.metadata.data.WriteCellData;
import com.alibaba.excel.metadata.property.ExcelContentProperty;
import com.mgi.pacs.common.core.exception.CloudPacsException;

//转化男女
public class GenderConverter implements Converter<Integer> {

    @Override
    public Integer convertToJavaData(ReadCellData<?> cellData, ExcelContentProperty contentProperty, GlobalConfiguration globalConfiguration) throws Exception {

        String value = cellData.getStringValue();
        value = StrUtil.trim(value);
        if ("男".equals(value)) {
            return 1;
        } else if ("女".equals(value)) {
            return 0;
        }
        throw new CloudPacsException("性别转换异常");
    }

    @Override
    public WriteCellData<?> convertToExcelData(Integer value, ExcelContentProperty contentProperty, GlobalConfiguration globalConfiguration) throws Exception {

        WriteCellData<String> cellData = new WriteCellData<>();
        if (value == 0) {
            cellData.setData("女");
        } else if (value == 1) {
            cellData.setData("男");
        }
        return cellData;
    }
}
