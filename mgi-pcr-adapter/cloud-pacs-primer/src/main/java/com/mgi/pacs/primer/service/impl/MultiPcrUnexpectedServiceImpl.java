package com.mgi.pacs.primer.service.impl;

import cn.hutool.core.collection.CollectionUtil;
import com.alibaba.excel.EasyExcel;
import com.alibaba.excel.read.listener.PageReadListener;
import com.alibaba.excel.read.metadata.ReadSheet;
import com.baomidou.mybatisplus.extension.service.impl.ServiceImpl;
import com.mgi.pacs.common.core.exception.CloudPacsException;
import com.mgi.pacs.primer.domain.pojo.MultiPcrEfficiency;
import com.mgi.pacs.primer.domain.pojo.MultiPcrUnexpected;
import com.mgi.pacs.primer.mapper.MultiPcrUnexpectedMapper;
import com.mgi.pacs.primer.service.IMultiPcrUnexpectedService;
import lombok.extern.slf4j.Slf4j;
import org.springframework.stereotype.Service;
import org.springframework.transaction.annotation.Transactional;
import org.springframework.web.multipart.MultipartFile;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

/**
 * <p>
 * 多重PCR（Polymerase Chain Reaction，聚合酶链反应）实验中出现的意外或非预期结果 服务实现类
 * </p>
 *
 * @author mabinbin
 * @since 2024-01-27
 */

@Slf4j
@Service
public class MultiPcrUnexpectedServiceImpl extends ServiceImpl<MultiPcrUnexpectedMapper, MultiPcrUnexpected> implements IMultiPcrUnexpectedService {

    /**
     * MultiPcr_unexpected
     * @param file 文件
     * @param experimentId 实验id
     */
    @Transactional(rollbackFor = Exception.class)
    @Override
    public void excelUpload(MultipartFile file, String experimentId) {

        log.info("开始解析 MultiPcr_unexpected {}", experimentId);
        // 创建一个列表来收集所有数据
        List<MultiPcrUnexpected> multiPcrUnexpecteds = new ArrayList<>();
        try {

            // 获取 Excel 文件中所有的sheet
            List<ReadSheet> sheets = EasyExcel.read(file.getInputStream())
                    .head(MultiPcrUnexpected.class)
                    .build()
                    .excelExecutor()
                    .sheetList();

            // 读取每个 sheet 并处理数据
            for (ReadSheet sheet : sheets) {
                String sheetName = sheet.getSheetName();
                EasyExcel.read(file.getInputStream(), MultiPcrUnexpected.class, new PageReadListener<MultiPcrUnexpected>(dataList -> {
                    for (MultiPcrUnexpected data : dataList) {
                        data.setExperiment_order(sheetName); // 设置 sheet 名称到 experiment_order 字段
                        multiPcrUnexpecteds.add(data); // 将处理后的数据添加到 experimentDataList
                    }
                })).sheet(sheet.getSheetNo()).doRead();
            }

            if (CollectionUtil.isEmpty(multiPcrUnexpecteds)) {
                throw new CloudPacsException("MultiPcr_unexpected 当前文件为空");
            }
        } catch (IOException e) {
            throw new CloudPacsException(e.getMessage());
        }
        log.info("解析 MultiPcr_unexpected {} 成功", experimentId);
        multiPcrUnexpecteds.forEach(e -> {
            e.setExperiment_id(experimentId);
        });

        this.saveBatch(multiPcrUnexpecteds);
    }
}
