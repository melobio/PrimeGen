package com.mgi.pacs.primer.controller;

import cn.hutool.core.lang.Assert;
import cn.hutool.core.util.ObjectUtil;
import com.baomidou.mybatisplus.core.conditions.query.LambdaQueryWrapper;
import com.mgi.pacs.common.core.exception.CloudPacsException;
import com.mgi.pacs.common.core.response.ServerResponseEntity;
import com.mgi.pacs.primer.domain.CommonBO;
import com.mgi.pacs.primer.domain.pojo.*;
import com.mgi.pacs.primer.mapper.CustomMapper;
import com.mgi.pacs.primer.service.*;
import io.swagger.annotations.Api;
import io.swagger.annotations.ApiOperation;
import io.swagger.annotations.ApiParam;
import lombok.extern.slf4j.Slf4j;
import org.apache.commons.collections4.map.LinkedMap;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.http.MediaType;
import org.springframework.transaction.annotation.Transactional;
import org.springframework.web.bind.annotation.*;
import org.springframework.web.multipart.MultipartFile;

import javax.validation.Valid;
import java.util.List;

/**
 * <p>
 * 引物表 前端控制器
 * </p>
 *
 * @author mabinbin
 * @since 2024-01-27
 */
@Api(tags = "primer相关接口")
@Slf4j
@RestController
@RequestMapping("/primer")
public class DtWqPrimerController {

    @Autowired
    private IDtWqPrimerService dtWqPrimerService;

    @Autowired
    private IMultiPcrSummaryService multiPcrSummaryService;

    @Autowired
    private IDimerInfoService dimerInfoService;

    @Autowired
    private IMultiPcrEfficiencyService multiPcrEfficiencyService;

    @Autowired
    private IMultiPcrUnexpectedService multiPcrUnexpectedService;

    @Autowired
    private CustomMapper customMapper;

    @Transactional(rollbackFor = Exception.class)
    @PostMapping(value = "/file/upload", consumes = MediaType.MULTIPART_FORM_DATA_VALUE)
    @ApiOperation(value = "上传primer数据excel", notes = "上传primer数据excel")
    public ServerResponseEntity<String> excelUpload(@RequestPart("files") MultipartFile[] files,
                                                  @ApiParam(name = "experimentId", value = "实验版本", example = "primer_DT1_1")
                                                  @RequestParam(value = "experimentId", required = false) String experimentId) {

        log.info("开始处理 {} 实验数据", experimentId);

        //TODO 重传现象待控制
        //没有文件上传
        if (files.length == 0) {
            return ServerResponseEntity.fail("请上传文件");
        }

        //需要检查文件列表中是否有文件版本，此时需要自定义文件版本
        String version = null;
        for (int i = 0; i < files.length; i++) {

            MultipartFile file = files[i];
            String filename = file.getOriginalFilename();
            Assert.isTrue(ObjectUtil.isNotEmpty(filename), () -> new CloudPacsException("filename is empty"));

            if (filename.contains("_primer_")) {
                version = dtWqPrimerService.getVersion(file);
            }
        }

        Assert.isTrue(version != null, () -> new CloudPacsException("need to upload primer table"));
        experimentId = version;

        //文件解析并保存数据库
        for (int i = 0; i < files.length; i++) {
            MultipartFile file = files[i];
            String filename = file.getOriginalFilename();

            log.info("正在处理文件解析文件：{}", filename);
            Assert.isTrue(ObjectUtil.isNotEmpty(filename), () -> new CloudPacsException("filename can not be empty"));

            //解析保存summary和depth
            if (filename.contains("MultiPcr_summary_depth")) {
                //避免重复上传
                int count = multiPcrSummaryService.count(new LambdaQueryWrapper<MultiPcrSummary>()
                        .eq(MultiPcrSummary::getExperiment_id, experimentId)
                        .last("for update"));
                if (count == 0) {
                    log.info("start save MultiPcr_summary_depth ...");
                    multiPcrSummaryService.excelUpload(file, experimentId);
                }
            }

            //解析保存dimer_info 有可能没有
            if (filename.contains("dimer_info")) {

                //get filename suffix
                String suffix = filename.substring(0, filename.indexOf("_dimer_info"));

                //避免重复上传
                int count = dimerInfoService.count(new LambdaQueryWrapper<DimerInfo>()
                        .eq(DimerInfo::getExperiment_id, experimentId)
                        .eq(DimerInfo::getExperiment_order, suffix)
                        .last("for update"));
                if (count == 0) {
                    dimerInfoService.excelUpload(file, experimentId, suffix);
                }
            }

            //解析MultiPcr_efficiency
            if (filename.contains("MultiPcr_efficiency")) {
                //避免重复上传
                int count = multiPcrEfficiencyService.count(new LambdaQueryWrapper<MultiPcrEfficiency>()
                        .eq(MultiPcrEfficiency::getExperiment_id, experimentId)
                        .last("for update"));
                if (count == 0) {
                    log.info("start save MultiPcr_efficiency_depth ...");
                    multiPcrEfficiencyService.excelUpload(file, experimentId);
                }
            }

            //解析MultiPcr_unexpected
            if (filename.contains("MultiPcr_unexpected")) {
                //避免重复上传
                int count = multiPcrUnexpectedService.count(new LambdaQueryWrapper<MultiPcrUnexpected>()
                        .eq(MultiPcrUnexpected::getExperiment_id, experimentId)
                        .last("for update"));
                if (count == 0) {
                    log.info("start save MultiPcr_unexpected_depth ...");
                    multiPcrUnexpectedService.excelUpload(file, experimentId);
                }
            }
        }

        return ServerResponseEntity.success(version);
    }

    @PostMapping(value = "/searchBySql")
    @ApiOperation(value = "根据sql查询结果", notes = "上传primer数据excel")
    public ServerResponseEntity<List<LinkedMap<String, Object>>> selectBySql(@Valid @RequestBody CommonBO commonBO) {

        List<LinkedMap<String, Object>> res = customMapper.executeSelect(commonBO.getSql());
        return ServerResponseEntity.success(res);
    }

}

