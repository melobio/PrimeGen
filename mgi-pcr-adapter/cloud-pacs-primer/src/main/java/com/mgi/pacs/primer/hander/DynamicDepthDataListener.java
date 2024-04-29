package com.mgi.pacs.primer.hander;

import com.alibaba.excel.context.AnalysisContext;
import com.alibaba.excel.event.AnalysisEventListener;
import com.mgi.pacs.primer.domain.pojo.MultiPcrDepth;

import java.math.BigDecimal;
import java.math.RoundingMode;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class DynamicDepthDataListener extends AnalysisEventListener<Map<Integer, String>> {

    private List<MultiPcrDepth> beans = new ArrayList<>();
    private List<String> headerNames = new ArrayList<>();

    @Override
    public void invokeHeadMap(Map<Integer, String> headMap, AnalysisContext context) {
        // 保存表头名称
        headMap.forEach((index, name) -> headerNames.add(name));
    }

    @Override
    public void invoke(Map<Integer, String> data, AnalysisContext context) {
        // 创建 DynamicBean 对象，并填充固定列
        String name = data.get(0); //当前行的name

        // 创建动态列的 Map
        data.forEach((index, value) -> {
            if (index > 3) { // 索引大于3的是动态列
                MultiPcrDepth bean = new MultiPcrDepth();
                bean.setName(name);
                bean.setExperiment_order(headerNames.get(index));
                bean.setDepth(new BigDecimal(value).setScale(2, RoundingMode.HALF_UP));
                // 将 bean 添加到列表中
                beans.add(bean);
            }
        });
    }

    @Override
    public void doAfterAllAnalysed(AnalysisContext context) {

    }

    public List<MultiPcrDepth> getBeans() {
        return beans;
    }
}
