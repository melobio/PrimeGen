package com.mgi.pacs.common.cache.constant;

/**
 * 缓存名字
 *
 * @author linzhunsheng
 * @date 2022/03/16
 */
public interface CacheNames {
    /**
     *
     * 参考CacheKeyPrefix
     * cacheNames 与 key 之间的默认连接字符
     */
    String UNION = "::";

    /**
     * key内部的连接字符（自定义）
     */
    String UNION_KEY = ":";
}
