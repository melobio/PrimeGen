package com.mgi.pacs.primer.service;

import com.baomidou.mybatisplus.extension.service.IService;
import com.mgi.pacs.primer.domain.pojo.DtWqPrimer;
import org.springframework.web.multipart.MultipartFile;

/**
 * <p>
 * 引物表 服务类
 * </p>
 *
 * @author mabinbin
 * @since 2024-01-27
 */
public interface IDtWqPrimerService extends IService<DtWqPrimer> {

    /**
     * 上传引物excel
     * @param file excel
     * @param experimentId 实验id
     */
    void excelUpload(MultipartFile file, String experimentId);

    /**
     * get version by primer table
     * @param file file
     * @return version
     */
    String getVersion(MultipartFile file);
}
