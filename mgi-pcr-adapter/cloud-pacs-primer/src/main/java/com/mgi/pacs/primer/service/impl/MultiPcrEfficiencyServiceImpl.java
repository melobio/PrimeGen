package com.mgi.pacs.primer.service.impl;

import cn.hutool.core.collection.CollectionUtil;
import com.alibaba.excel.EasyExcel;
import com.alibaba.excel.read.listener.PageReadListener;
import com.alibaba.excel.read.metadata.ReadSheet;
import com.baomidou.mybatisplus.extension.service.impl.ServiceImpl;
import com.mgi.pacs.common.core.exception.CloudPacsException;
import com.mgi.pacs.primer.domain.pojo.MultiPcrEfficiency;
import com.mgi.pacs.primer.mapper.MultiPcrEfficiencyMapper;
import com.mgi.pacs.primer.service.IMultiPcrEfficiencyService;
import lombok.extern.slf4j.Slf4j;
import org.springframework.stereotype.Service;
import org.springframework.transaction.annotation.Transactional;
import org.springframework.web.multipart.MultipartFile;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

/**
 * <p>
 * 多重PCR(聚合酶链反应)实验中的效率 服务实现类
 * </p>
 *
 * @author mabinbin
 * @since 2024-01-27
 */
@Slf4j
@Service
public class MultiPcrEfficiencyServiceImpl extends ServiceImpl<MultiPcrEfficiencyMapper, MultiPcrEfficiency> implements IMultiPcrEfficiencyService {

    /**
     * 解析MultiPcr_efficiency
     * @param file 文件
     * @param experimentId 实验id
     */
    @Transactional(rollbackFor = Exception.class)
    @Override
    public void excelUpload(MultipartFile file, String experimentId) {


        log.info("开始解析 MultiPcr_efficiency {}", experimentId);
        // 创建一个列表来收集所有数据
        List<MultiPcrEfficiency> experimentDataList = new ArrayList<>();
        try {

            // 获取 Excel 文件中所有的sheet
            List<ReadSheet> sheets = EasyExcel.read(file.getInputStream())
                    .head(MultiPcrEfficiency.class)
                    .build()
                    .excelExecutor()
                    .sheetList();

            // 读取每个 sheet 并处理数据
            for (ReadSheet sheet : sheets) {
                String sheetName = sheet.getSheetName();
                EasyExcel.read(file.getInputStream(), MultiPcrEfficiency.class, new PageReadListener<MultiPcrEfficiency>(dataList -> {
                    for (MultiPcrEfficiency data : dataList) {
                        data.setExperiment_order(sheetName); // 设置 sheet 名称到 experiment_order 字段
                        experimentDataList.add(data); // 将处理后的数据添加到 experimentDataList
                    }
                })).sheet(sheet.getSheetNo()).doRead();
            }

            if (CollectionUtil.isEmpty(experimentDataList)) {
                throw new CloudPacsException("MultiPcr_efficiency 当前文件为空");
            }
        } catch (IOException e) {
            throw new CloudPacsException(e.getMessage());
        }
        log.info("解析 MultiPcr_efficiency {} 成功", experimentId);
        experimentDataList.forEach(e -> {
            e.setExperiment_id(experimentId);
        });

        this.saveBatch(experimentDataList);
    }
}
