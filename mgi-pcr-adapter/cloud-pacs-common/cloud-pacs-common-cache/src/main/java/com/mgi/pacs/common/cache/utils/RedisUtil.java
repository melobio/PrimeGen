package com.mgi.pacs.common.cache.utils;

import cn.hutool.core.collection.CollUtil;
import cn.hutool.core.util.StrUtil;
import com.mgi.pacs.common.cache.constant.CacheNames;
import com.mgi.pacs.common.core.exception.CloudPacsException;
import com.mgi.pacs.common.core.response.ResponseEnum;
import lombok.extern.slf4j.Slf4j;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.data.redis.core.BoundSetOperations;
import org.springframework.data.redis.core.HashOperations;
import org.springframework.data.redis.core.RedisTemplate;
import org.springframework.data.redis.core.StringRedisTemplate;
import org.springframework.data.redis.core.script.DefaultRedisScript;
import org.springframework.stereotype.Component;

import javax.annotation.PostConstruct;
import java.util.*;
import java.util.concurrent.TimeUnit;
import java.util.stream.Collectors;

/**
 * redis缓存工具类
 *
 * @author linzhunsheng
 * @date 2022/03/16
 */
@Slf4j
@Component
public class RedisUtil {

    private static RedisUtil redisUtil;

    @Autowired
    private RedisTemplate redisTemplate;

    @Autowired
    private StringRedisTemplate stringRedisTemplate;

    @PostConstruct
    public void init() {
        redisUtil = this;
        redisUtil.redisTemplate = this.redisTemplate;
        redisUtil.stringRedisTemplate = this.stringRedisTemplate;
    }

    // =============================common============================
    /**
     * 指定缓存失效时间
     * @param key 键
     * @param time 时间(秒)
     * @return 是否成功
     */
    public static Boolean expire(String key, long time) {
        return expire(key, time, null);
    }

    /**
     * 设置有效时间
     *
     * @param key Redis键
     * @param timeout 超时时间
     * @param unit 时间单位
     * @return true=设置成功；false=设置失败
     */
    public static Boolean expire(String key, long timeout, TimeUnit unit) {
        if (key.contains(StrUtil.SPACE)) {
            throw new CloudPacsException(ResponseEnum.EXCEPTION);
        }
        try {
            if (timeout > 0 && unit != null) {
                redisUtil.redisTemplate.expire(key, timeout, unit);
            } else if (timeout > 0) {
                redisUtil.redisTemplate.expire(key, timeout, TimeUnit.SECONDS);
            }
            return Boolean.TRUE;
        } catch (Exception e) {
            log.error("Set expire error: {}", e.getMessage());
            return Boolean.FALSE;
        }
    }

    /**
     * 根据key 获取过期时间
     * @param key 键 不能为null
     * @return 时间(秒) 返回-1代表为永久有效 失效时间为0，说明该主键未设置失效时间（失效时间默认为-1）
     */
    public static Long getExpire(String key) {
        if (key.contains(StrUtil.SPACE)) {
            throw new CloudPacsException(ResponseEnum.EXCEPTION);
        }
        return redisUtil.redisTemplate.getExpire(key, TimeUnit.SECONDS);
    }

    /**
     * 判断key是否存在
     * @param key 键
     * @return true 存在 false 不存在
     */
    public static Boolean hasKey(String key) {
        if (key.contains(StrUtil.SPACE)) {
            throw new CloudPacsException(ResponseEnum.EXCEPTION);
        }
        try {
            return redisUtil.redisTemplate.hasKey(key);
        }
        catch (Exception e) {
            log.error("Error getting hasKey: {}", e.getMessage());
            return Boolean.FALSE;
        }
    }

    /**
     * 删除缓存
     * @param key 可以传一个值 或多个
     */
    @SuppressWarnings("unchecked")
    public static void del(String... key) {
        if (key != null && key.length > 0) {
            for (String s : key) {
                if (s.contains(StrUtil.SPACE)) {
                    throw new CloudPacsException(ResponseEnum.EXCEPTION);
                }
            }

            if (key.length == 1) {
                redisUtil.redisTemplate.delete(key[0]);
            }
            else {
                redisUtil.redisTemplate.delete(Arrays.asList(key));
            }
        }
    }

    // ============================String, Integer, Object=============================
    /**
     * 普通缓存获取
     * @param key 键
     * @return 值
     */
    @SuppressWarnings("unchecked")
    public static <T> T getCacheObject(String key) {
        if (key.contains(StrUtil.SPACE)) {
            throw new CloudPacsException(ResponseEnum.EXCEPTION);
        }
        return (T) redisUtil.redisTemplate.opsForValue().get(key);
    }

    /**
     * 缓存基本的对象，Integer、String、实体类等
     *
     * @param key 缓存的键值
     * @param value 缓存的值
     */
    public static <T> boolean setCacheObject(String key, T value) {
        return setCacheObject(key, value, 0);
    }

    /**
     * 普通缓存放入并设置时间
     * @param key 键
     * @param value 值
     * @param time 时间(秒) time要大于0 如果time小于等于0 将设置无限期
     * @return true成功 false 失败
     */
    public static <T> boolean setCacheObject(String key, T value, long time) {
        return setCacheObject(key, value, time, null);
    }

    /**
     * 缓存基本的对象，Integer、String、实体类等
     *
     * @param key 缓存的键值
     * @param value 缓存的值
     * @param timeout 时间
     * @param timeUnit 时间颗粒度
     */
    public static <T> boolean setCacheObject(String key, T value, long timeout, TimeUnit timeUnit) {
        if (key.contains(StrUtil.SPACE)) {
            throw new CloudPacsException(ResponseEnum.EXCEPTION);
        }
        try {
            if (timeout > 0 && timeUnit != null) {
                redisUtil.redisTemplate.opsForValue().set(key, value, timeout, timeUnit);
            } else if (timeout > 0) {
                redisUtil.redisTemplate.opsForValue().set(key, value, timeout, TimeUnit.SECONDS);
            } else {
                redisUtil.redisTemplate.opsForValue().set(key, value);
            }
            return true;
        }
        catch (Exception e) {
            log.error("Redis opsForValue error: {}", e.getMessage());
            return false;
        }
    }

    /**
     * 递增 此时value值必须为int类型 否则报错
     * @param key 键
     * @param delta 要增加几(大于0)
     * @return 自增后的值
     */
    public static Long incr(String key, long delta) {
        if (key.contains(StrUtil.SPACE)) {
            throw new CloudPacsException(ResponseEnum.EXCEPTION);
        }
        if (delta < 0) {
            throw new RuntimeException("递增因子必须大于0");
        }
        return redisUtil.stringRedisTemplate.opsForValue().increment(key, delta);
    }

    /**
     * 递减
     * @param key 键
     * @param delta 要减少几(小于0)
     * @return 自减后的值
     */
    public static Long decr(String key, long delta) {
        if (key.contains(StrUtil.SPACE)) {
            throw new CloudPacsException(ResponseEnum.EXCEPTION);
        }
        if (delta < 0) {
            throw new RuntimeException("递减因子必须小于0");
        }
        return redisUtil.stringRedisTemplate.opsForValue().increment(key, -delta);
    }

    public static boolean setLongValue(String key, Long value, long time) {
        if (key.contains(StrUtil.SPACE)) {
            throw new CloudPacsException(ResponseEnum.EXCEPTION);
        }
        try {
            if (time > 0) {
                redisUtil.stringRedisTemplate.opsForValue().set(key, String.valueOf(value), time, TimeUnit.SECONDS);
            }
            else {
                redisUtil.stringRedisTemplate.opsForValue().set(key, String.valueOf(value));
            }
            return true;
        }
        catch (Exception e) {
            log.error("setLongValue() error: {}", e.getMessage());
            return false;
        }
    }

    /**
     * 普通缓存获取
     * @param key 键
     * @return 值
     */
    public static Long getLongValue(String key) {
        if (key == null) {
            return null;
        }
        if (key.contains(StrUtil.SPACE)) {
            throw new CloudPacsException(ResponseEnum.EXCEPTION);
        }
        String result = redisUtil.stringRedisTemplate.opsForValue().get(key);
        if (result == null) {
            return null;
        }
        return Long.valueOf(result);
    }

    /**
     * 批量删除缓存
     * @param keys
     */
    public static void deleteBatch(List<String> keys) {
        if (CollUtil.isEmpty(keys)) {
            return;
        }
        for (String key : keys) {
            if (key.contains(StrUtil.SPACE)) {
                throw new CloudPacsException(ResponseEnum.EXCEPTION);
            }
        }
        redisUtil.redisTemplate.delete(keys);
    }

    /**
     * 批量删除缓存
     * @param cacheName 缓存名
     * @param cacheKeys 缓存key
     */
    public static void deleteBatch(String cacheName, List<?> cacheKeys) {
        if (StrUtil.isBlank(cacheName) || CollUtil.isEmpty(cacheKeys)) {
            return;
        }
        List<String> strCacheKeys = cacheKeys.stream().map(String::valueOf).collect(Collectors.toList());
        List<String> keys = new ArrayList<>();
        for (String cacheKey : strCacheKeys) {
            String key = cacheName + CacheNames.UNION + cacheKey;
            keys.add(key);
            if (key.contains(StrUtil.SPACE)) {
                throw new CloudPacsException(ResponseEnum.EXCEPTION);
            }
        }
        redisUtil.redisTemplate.delete(keys);
    }

    /**
     * 比较和删除标记，原子性
     * @return 是否成功
     */
    public static boolean cad(String key, String value) {

        if (key.contains(StrUtil.SPACE) || value.contains(StrUtil.SPACE)) {
            throw new CloudPacsException(ResponseEnum.EXCEPTION);
        }

        String script = "if redis.call('get', KEYS[1]) == ARGV[1] then return redis.call('del', KEYS[1]) else return 0 end";

        //通过lure脚本原子验证令牌和删除令牌
        Long result = redisUtil.stringRedisTemplate.execute(new DefaultRedisScript<Long>(script, Long.class),
                Collections.singletonList(key),
                value);

        return !Objects.equals(result, 0L);
    }

    // ============================List=============================
    /**
     * 缓存List数据
     *
     * @param key      缓存的键值
     * @param dataList 待缓存的List数据
     * @return 缓存的对象
     */
    public static <T> long setCacheList(String key, List<T> dataList) {
        Long count = redisUtil.redisTemplate.opsForList().rightPushAll(key, dataList);
        return count == null ? 0 : count;
    }

    /**
     * 获得缓存的list对象
     *
     * @param key 缓存的键值
     * @return 缓存键值对应的数据
     */
    public static List<Object> getCacheList(String key) {
        return redisUtil.redisTemplate.opsForList().range(key, 0, -1);
    }

    // ============================Set=============================
    /**
     * 缓存Set
     *
     * @param key 缓存键值
     * @param dataSet 缓存的数据
     * @return 缓存数据的对象
     */
    public static <T> BoundSetOperations<String, Object> setCacheSet(String key, Set<T> dataSet) {
        BoundSetOperations<String, Object> setOperation = redisUtil.redisTemplate.boundSetOps(key);
        Iterator<T> it = dataSet.iterator();
        while (it.hasNext()) {
            setOperation.add(it.next());
        }
        return setOperation;
    }

    /**
     * 缓存Set
     *
     * @param key 缓存键值
     * @param data 缓存的数据
     * @return 缓存数据的对象
     */
    public static <T> BoundSetOperations<String, Object> setCacheSet(String key, T data) {
        BoundSetOperations<String, Object> setOperation = redisUtil.redisTemplate.boundSetOps(key);
        setOperation.add(data);
        return setOperation;
    }

    /**
     * 获得缓存的set
     *
     * @param key
     * @return
     */
    public static <T> Set<Object> getCacheSet(String key) {
        return redisUtil.redisTemplate.opsForSet().members(key);
    }

    /**
     * 是否包含val
     *
     * @param key
     * @param data
     * @return
     */
    public static boolean cacheSetContains(String key, Object data) {
        return redisUtil.redisTemplate.opsForSet().isMember(key, data);
    }

    /**
     * 移除值为value的
     *
     * @param key    键
     * @param values 值 可以是多个
     * @return 移除的个数
     */
    public static long cacheSetRemove(String key, Object... values) {
        try {
            Long count = redisUtil.redisTemplate.opsForSet().remove(key, values);
            return count;
        } catch (Exception e) {
            e.printStackTrace();
            return 0;
        }
    }

    // ============================Map=============================
    /**
     * 缓存Map
     *
     * @param key
     * @param dataMap
     */
    public static <T> void setCacheMap(String key, Map<String, T> dataMap) {
        if (dataMap != null) {
            redisUtil.redisTemplate.opsForHash().putAll(key, dataMap);
        }
    }

    /**
     * 获得缓存的Map
     *
     * @param key
     * @return
     */
    public static <T> Map<Object, Object> getCacheMap(String key) {
        return redisUtil.redisTemplate.opsForHash().entries(key);
    }

    /**
     * 指定键，删除对应的值
     *
     * @param key Redis键
     * @param hKey Hash键
     * @param <T> 值
     */
    public static <T> void delCacheMap(String key, String hKey) {
        redisUtil.redisTemplate.opsForHash().delete(key, hKey);
    }

    /**
     * 往Hash中存入数据
     *
     * @param key Redis键
     * @param hKey Hash键
     * @param value 值
     */
    public static <T> void setCacheMapValue(final String key, final String hKey, final T value) {
        redisUtil.redisTemplate.opsForHash().put(key, hKey, value);
    }

    /**
     * 获取Hash中的数据
     *
     * @param key Redis键
     * @param hKey Hash键
     * @return Hash中的对象
     */
    public static <T> T getCacheMapValue(final String key, final String hKey) {
        HashOperations<String, String, T> opsForHash = redisUtil.redisTemplate.opsForHash();
        return opsForHash.get(key, hKey);
    }

    /**
     * 获取多个Hash中的数据
     *
     * @param key Redis键
     * @param hKeys Hash键集合
     * @return Hash对象集合
     */
    public static <T> List<Object> getMultiCacheMapValue(final String key, final Collection<Object> hKeys) {
        return redisUtil.redisTemplate.opsForHash().multiGet(key, hKeys);
    }

    /**
     * 获得缓存的基本对象列表
     *
     * @param pattern 字符串前缀
     * @return 对象列表
     */
    public static Collection<String> keys(final String pattern) {
        return redisUtil.redisTemplate.keys(pattern);
    }
}
