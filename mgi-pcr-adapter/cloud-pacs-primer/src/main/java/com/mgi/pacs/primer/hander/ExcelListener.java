package com.mgi.pacs.primer.hander;

import com.alibaba.excel.context.AnalysisContext;
import com.alibaba.excel.read.listener.ReadListener;
import com.mgi.pacs.primer.domain.BgiFamilyOrderReq;
import lombok.Data;
import lombok.extern.slf4j.Slf4j;

import java.util.ArrayList;
import java.util.List;

@Slf4j
@Data
public class ExcelListener implements ReadListener {

    private List dataList = new ArrayList<>();

    @Override
    public void invoke(Object data, AnalysisContext context) {
        dataList.add(data);
    }

    @Override
    public void doAfterAllAnalysed(AnalysisContext context) {

    }

    public List getDataList() {
        return dataList;
    }
}