package com.mgi.pacs.primer.service.impl;

import cn.hutool.core.collection.CollectionUtil;
import com.alibaba.excel.annotation.ExcelProperty;
import com.baomidou.mybatisplus.extension.service.impl.ServiceImpl;
import com.mgi.pacs.primer.domain.pojo.DimerInfo;
import com.mgi.pacs.primer.mapper.DimerInfoMapper;
import com.mgi.pacs.primer.service.IDimerInfoService;
import com.mgi.pacs.primer.utils.MgiExcelUtil;
import lombok.extern.slf4j.Slf4j;
import org.springframework.stereotype.Service;
import org.springframework.transaction.annotation.Transactional;
import org.springframework.web.multipart.MultipartFile;

import java.util.List;

/**
 * <p>
 * 二聚体信息 服务实现类
 * </p>
 *
 * @author mabinbin
 * @since 2024-01-27
 */
@Slf4j
@Service
public class DimerInfoServiceImpl extends ServiceImpl<DimerInfoMapper, DimerInfo> implements IDimerInfoService {

    /**
     * 解析保存dimer_info 有可能没有
     * @param file 文件
     * @param experimentId 实验id
     * @Param suffix experiment order
     */
    @Transactional(rollbackFor = Exception.class)
    @Override
    public void excelUpload(MultipartFile file, String experimentId, String suffix) {

        log.info("正在解析DimerInfo {} ...", experimentId);
        List<DimerInfo> dimerInfos = MgiExcelUtil.getListBySheetName(file, null, DimerInfo.class);
        log.info("解析DimerInfo成功 {} ...", experimentId);
        if (CollectionUtil.isEmpty(dimerInfos)) {
            return;
        }
        dimerInfos.forEach(d -> {
            d.setExperiment_id(experimentId);
            d.setExperiment_order(suffix);
        });
        this.saveBatch(dimerInfos);
    }
}
