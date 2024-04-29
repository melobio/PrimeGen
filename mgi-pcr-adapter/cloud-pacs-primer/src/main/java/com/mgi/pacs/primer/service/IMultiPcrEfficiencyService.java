package com.mgi.pacs.primer.service;

import com.baomidou.mybatisplus.extension.service.IService;
import com.mgi.pacs.primer.domain.pojo.MultiPcrEfficiency;
import org.springframework.web.multipart.MultipartFile;

/**
 * <p>
 * 多重PCR(聚合酶链反应)实验中的效率 服务类
 * </p>
 *
 * @author mabinbin
 * @since 2024-01-27
 */
public interface IMultiPcrEfficiencyService extends IService<MultiPcrEfficiency> {

    /**
     * 解析MultiPcr_efficiency
     * @param file 文件
     * @param experimentId 实验id
     */
    void excelUpload(MultipartFile file, String experimentId);
}
