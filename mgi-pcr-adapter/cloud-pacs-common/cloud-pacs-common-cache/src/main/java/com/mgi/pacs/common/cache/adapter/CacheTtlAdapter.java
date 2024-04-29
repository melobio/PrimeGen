package com.mgi.pacs.common.cache.adapter;

import com.mgi.pacs.common.cache.bo.CacheNameWithTtlBO;

import java.util.List;

/**
 * 实现该接口之后，根据缓存的cacheName和ttl将缓存进行过期
 *
 * @author linzhunsheng
 * @date 2022/03/16
 */
public interface CacheTtlAdapter {

    /**
     * 根据缓存的cacheName和ttl将缓存进行过期
     * @return 需要独立设置过期时间的缓存列表
     */
    List<CacheNameWithTtlBO> listCacheNameWithTtl();
}
