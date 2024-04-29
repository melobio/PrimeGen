package com.mgi.pacs.primer.service;

import com.baomidou.mybatisplus.extension.service.IService;
import com.mgi.pacs.primer.domain.pojo.MultiPcrSummary;
import org.springframework.web.multipart.MultipartFile;

/**
 * <p>
 * 多重PCR统计数据 服务类
 * </p>
 *
 * @author mabinbin
 * @since 2024-01-27
 */
public interface IMultiPcrSummaryService extends IService<MultiPcrSummary> {

    /**
     * 解析保存summary和depth
     * @param file 文件
     * @param experimentId 实验id
     */
    void excelUpload(MultipartFile file, String experimentId);
}
