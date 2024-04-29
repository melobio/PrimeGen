package com.mgi.pacs.primer.service;

import com.baomidou.mybatisplus.extension.service.IService;
import com.mgi.pacs.primer.domain.pojo.DimerInfo;
import org.springframework.web.multipart.MultipartFile;

/**
 * <p>
 * 二聚体信息 服务类
 * </p>
 *
 * @author mabinbin
 * @since 2024-01-27
 */
public interface IDimerInfoService extends IService<DimerInfo> {

    /**
     * 解析保存dimer_info 有可能没有
     * @param file 文件
     * @param experimentId 实验id
     * @Param suffix experiment order
     */
    void excelUpload(MultipartFile file, String experimentId, String suffix);
}
