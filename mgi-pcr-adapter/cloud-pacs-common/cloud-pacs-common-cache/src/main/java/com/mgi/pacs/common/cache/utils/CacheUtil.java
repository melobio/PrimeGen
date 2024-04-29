package com.mgi.pacs.common.cache.utils;


import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;

/**
 * 自定义本地缓存工具类
 *
 * @author linzhunsheng
 * @date 2022/03/16
 */
public class CacheUtil {

    /**
     * 本地缓存保存最大的数量
     */
    private static final int MAX_MAP_SIZE = 200;

    private static Map<String, CacheEntity> cache = new ConcurrentHashMap<>();

    public static void put(String key, Object value, long duration) {
        capacityCheck();
        cache.put(key, new CacheEntity(value, duration));
    }

    public static Object get(String key) {
        CacheEntity o = cache.get(key);
        if (o == null) {
            return null;
        }
        if (o.getExpirationTime() != 0) {
            if (o.getExpirationTime() < System.currentTimeMillis()) {
                return null;
            }
        }
        return o.getData();
    }

    private static void capacityCheck() {
        int size = cache.size();
        if (size > MAX_MAP_SIZE) {
            cache.clear();
        }
    }

    private static class CacheEntity {
        private long expirationTime;
        private Object data;

        public CacheEntity(Object data, long duration) {
            this.data = data;
            if (duration <= 0) {
                this.expirationTime = 0;
            } else {
                this.expirationTime = System.currentTimeMillis() + duration;
            }
        }

        public long getExpirationTime() {
            return expirationTime;
        }

        public Object getData() {
            return data;
        }
    }
}
