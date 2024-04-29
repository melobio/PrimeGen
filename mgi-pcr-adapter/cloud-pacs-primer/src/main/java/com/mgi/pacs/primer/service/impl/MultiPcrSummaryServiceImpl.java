package com.mgi.pacs.primer.service.impl;

import com.alibaba.excel.EasyExcel;
import com.baomidou.mybatisplus.extension.service.impl.ServiceImpl;
import com.mgi.pacs.common.core.exception.CloudPacsException;
import com.mgi.pacs.primer.domain.pojo.MultiPcrDepth;
import com.mgi.pacs.primer.domain.pojo.MultiPcrSummary;
import com.mgi.pacs.primer.hander.DynamicDepthDataListener;
import com.mgi.pacs.primer.hander.ExcelListener;
import com.mgi.pacs.primer.mapper.MultiPcrSummaryMapper;
import com.mgi.pacs.primer.service.IMultiPcrDepthService;
import com.mgi.pacs.primer.service.IMultiPcrSummaryService;
import com.mgi.pacs.primer.utils.MgiExcelUtil;
import lombok.extern.slf4j.Slf4j;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Service;
import org.springframework.transaction.annotation.Transactional;
import org.springframework.web.multipart.MultipartFile;

import java.io.IOException;
import java.util.List;

/**
 * <p>
 * 多重PCR统计数据 服务实现类
 * </p>
 *
 * @author mabinbin
 * @since 2024-01-27
 */

@Slf4j
@Service
public class MultiPcrSummaryServiceImpl extends ServiceImpl<MultiPcrSummaryMapper, MultiPcrSummary> implements IMultiPcrSummaryService {

    @Autowired
    private IMultiPcrDepthService multiPcrDepthService;

    /**
     * 解析保存summary和depth
     * @param file 文件
     * @param experimentId 实验id
     */
    @Transactional(rollbackFor = Exception.class)
    @Override
    public void excelUpload(MultipartFile file, String experimentId) {

        log.info("开始保存summary数据 {} ...", experimentId);
        List<MultiPcrSummary> summaryList = MgiExcelUtil.getListBySheetName(file, "summary", MultiPcrSummary.class);
        summaryList.forEach((s) -> {
            s.setExperiment_id(experimentId);
        });
        log.info("开始保存summary数据成功 {} ...", experimentId);

        log.info("开始保存depth数据 {} ...", experimentId);
        List<MultiPcrDepth> depths = null;
        try {
            DynamicDepthDataListener dynamicDepthDataListener = new DynamicDepthDataListener();
            EasyExcel.read(file.getInputStream(), dynamicDepthDataListener).sheet("amp_unique_region_depth").doRead();
            depths = dynamicDepthDataListener.getBeans();
            if (depths.size() == 0) {
                throw new CloudPacsException("请检查excel中是否写入数据");
            }
        } catch (IOException e) {
            throw new CloudPacsException("解析excel失败");
        }

        depths.forEach(d -> {
            d.setExperiment_id(experimentId);
        });
        log.info("开始保存depth数据成功 {} ...", experimentId);

        this.saveBatch(summaryList);
        multiPcrDepthService.saveBatch(depths);
    }


}
