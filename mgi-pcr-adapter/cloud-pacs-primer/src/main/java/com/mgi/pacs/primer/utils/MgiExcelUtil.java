package com.mgi.pacs.primer.utils;

import com.alibaba.excel.EasyExcel;
import com.mgi.pacs.common.core.exception.CloudPacsException;
import com.mgi.pacs.primer.domain.pojo.MultiPcrSummary;
import com.mgi.pacs.primer.hander.ExcelListener;
import org.springframework.web.multipart.MultipartFile;

import java.io.IOException;
import java.util.List;

public class MgiExcelUtil {

    /**
     * 获取excel对应sheet的列表
     * @param file 文件
     * @param sheetName sheetName
     * @return 列表
     */
    public static List getListBySheetName(MultipartFile file, String sheetName, Class clazz) {
        List excelResDataList = null;

        try {
            ExcelListener excelListener = new ExcelListener();

            if (sheetName != null) {
                EasyExcel.read(file.getInputStream(), clazz, excelListener).sheet(sheetName).doRead();
            } else {
                EasyExcel.read(file.getInputStream(), clazz, excelListener).sheet().doRead();
            }
            excelResDataList = excelListener.getDataList();

            if (excelResDataList.size() == 0) {
                throw new CloudPacsException("请检查excel中是否写入数据");
            }
        } catch (IOException e) {
            throw new CloudPacsException("解析excel失败");
        }

        return excelResDataList;
    }
}
