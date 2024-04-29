package com.mgi.pacs.common.core.utils;

import cn.hutool.core.util.ObjectUtil;
import cn.hutool.core.util.StrUtil;
import com.github.pagehelper.PageHelper;
import com.github.pagehelper.PageInfo;
import com.mgi.pacs.common.core.exception.CloudPacsException;
import com.mgi.pacs.common.core.page.PageDomain;
import com.mgi.pacs.common.core.page.PageVO;
import com.mgi.pacs.common.core.page.TableSupport;

import java.util.List;

/**
 * 分页工具类
 *
 * @author linzhunsheng
 * @date 2022/12/13
 */
public class PageUtils extends PageHelper {

    /**
     * 设置请求分页数据
     */
    public static void startPage() {
        PageDomain pageDomain = TableSupport.buildPageRequest();
        Integer pageNum = pageDomain.getPageNum();
        Integer pageSize = pageDomain.getPageSize();
        if (ObjectUtil.isNotNull(pageNum) && ObjectUtil.isNotNull(pageSize)) {
            String orderField = pageDomain.getOrderBy();
            // 仅支持字母、数字、下划线、空格、逗号、小数点（支持多个字段排序）
            String SQL_PATTERN = "[a-zA-Z0-9_\\ \\,\\.]+";
            if (StrUtil.isNotEmpty(orderField) && !orderField.matches(SQL_PATTERN)) {
                throw new CloudPacsException("参数不符合规范，不能进行查询");
            }

            PageHelper.startPage(pageNum, pageSize, orderField);
        }
    }

    /**
     * 响应请求分页数据
     */
    @SuppressWarnings({"rawtypes", "unchecked"})
    public static <T> PageVO<T> getPageVO(List<?> list) {
        PageVO pageVO = new PageVO();
        pageVO.setRows(list);
        pageVO.setTotal(new PageInfo(list).getTotal());
        return pageVO;
    }
}


