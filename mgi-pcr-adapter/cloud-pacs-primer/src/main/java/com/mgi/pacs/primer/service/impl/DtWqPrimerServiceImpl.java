package com.mgi.pacs.primer.service.impl;

import cn.hutool.core.util.ObjectUtil;
import com.alibaba.excel.EasyExcel;
import com.baomidou.mybatisplus.core.conditions.query.LambdaQueryWrapper;
import com.baomidou.mybatisplus.extension.service.impl.ServiceImpl;
import com.mgi.pacs.common.core.exception.CloudPacsException;
import com.mgi.pacs.primer.domain.pojo.DtWqPrimer;
import com.mgi.pacs.primer.hander.DtWqPrimerExcelListener;
import com.mgi.pacs.primer.mapper.DtWqPrimerMapper;
import com.mgi.pacs.primer.service.IDtWqPrimerService;
import lombok.extern.slf4j.Slf4j;
import org.springframework.stereotype.Service;
import org.springframework.transaction.annotation.Transactional;
import org.springframework.web.multipart.MultipartFile;

import java.io.IOException;
import java.util.List;

/**
 * <p>
 * 引物表 服务实现类
 * </p>
 *
 * @author mabinbin
 * @since 2024-01-27
 */
@Slf4j
@Service
public class DtWqPrimerServiceImpl extends ServiceImpl<DtWqPrimerMapper, DtWqPrimer> implements IDtWqPrimerService {

    /**
     * 上传引物excel
     * @param file excel
     * @param experimentId 实验id
     */
    @Transactional(rollbackFor = Exception.class)
    @Override
    public void excelUpload(MultipartFile file, String experimentId) {

        log.info("开始解析{} DtWqPrimer", experimentId);
        List<DtWqPrimer> dtWqPrimers = null;

        try {
            DtWqPrimerExcelListener dtWqPrimerExcelListener = new DtWqPrimerExcelListener();
            EasyExcel.read(file.getInputStream(), DtWqPrimer.class, dtWqPrimerExcelListener).sheet().doRead();
            dtWqPrimers = dtWqPrimerExcelListener.getDataList();

            if (dtWqPrimers.size() == 0) {
                throw new CloudPacsException("请检查excel中是否写入数据");
            }
        } catch (IOException e) {
            throw new CloudPacsException("解析excel失败");
        }

        //保存列表数据
        dtWqPrimers.forEach((dt) -> {
            dt.setExperiment_id(experimentId);
        });

        //保存数据库
        this.saveBatch(dtWqPrimers);
        log.info("解析{} DtWqPrimer成功", experimentId);
    }

    /**
     * get version by primer table
     * @param file file
     * @return version
     */
    @Transactional(rollbackFor = Exception.class)
    @Override
    public String getVersion(MultipartFile file) {

        log.info("start parse DtWqPrimer get version");
        String experimentId = "primer-";
        List<DtWqPrimer> dtWqPrimers = null;

        try {
            DtWqPrimerExcelListener dtWqPrimerExcelListener = new DtWqPrimerExcelListener();
            EasyExcel.read(file.getInputStream(), DtWqPrimer.class, dtWqPrimerExcelListener).sheet().doRead();
            dtWqPrimers = dtWqPrimerExcelListener.getDataList();

            if (dtWqPrimers.size() == 0) {
                throw new CloudPacsException("please check primer table have data ?");
            }
        } catch (IOException e) {
            throw new CloudPacsException("parse excel fail");
        }

        DtWqPrimer dtWqPrimer = dtWqPrimers.get(0);
        String version = ObjectUtil.isEmpty(dtWqPrimer.getVersion()) ? "v1" : dtWqPrimer.getVersion();
        experimentId = experimentId + dtWqPrimer.getChr_name1() + "-" + version;

        int count = this.count(new LambdaQueryWrapper<DtWqPrimer>()
                .eq(DtWqPrimer::getExperiment_id, experimentId)
                .last("for update"));
        if (count == 0) {
            //保存列表数据
            String finalVersion = experimentId;
            dtWqPrimers.forEach((dt) -> {
                dt.setExperiment_id(finalVersion);
            });

            //保存数据库
            this.saveBatch(dtWqPrimers);
        }

        return experimentId;
    }
}
