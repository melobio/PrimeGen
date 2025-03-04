package com.mgi.pacs.primer.service;

import com.baomidou.mybatisplus.extension.service.IService;
import com.mgi.pacs.primer.domain.pojo.MultiPcrUnexpected;
import org.springframework.web.multipart.MultipartFile;

/**
 * <p>
 * 多重PCR（Polymerase Chain Reaction，聚合酶链反应）实验中出现的意外或非预期结果 服务类
 * </p>
 *
 * @author mabinbin
 * @since 2024-01-27
 */
public interface IMultiPcrUnexpectedService extends IService<MultiPcrUnexpected> {

    /**
     * MultiPcr_unexpected
     * @param file 文件
     * @param experimentId 实验id
     */
    void excelUpload(MultipartFile file, String experimentId);
}
