package com.mgi.pacs.common.core.utils;


import cn.hutool.core.codec.Base64;
import cn.hutool.core.util.StrUtil;
import cn.hutool.json.JSONObject;
import cn.hutool.json.JSONUtil;
import io.jsonwebtoken.Claims;
import io.jsonwebtoken.Jws;
import io.jsonwebtoken.Jwts;
import lombok.extern.slf4j.Slf4j;
import org.apache.shiro.SecurityUtils;
import org.springframework.core.io.ClassPathResource;
import org.springframework.core.io.Resource;
import org.springframework.util.FileCopyUtils;
import sun.misc.BASE64Decoder;

import java.security.KeyFactory;
import java.security.PublicKey;
import java.security.spec.X509EncodedKeySpec;
import java.util.Optional;

/**
 * Jwt工具类
 *
 * @author linzhunsheng
 * @date 2022/03/16
 */
@Slf4j
public class JwtUtils {

    public static final long expire = 30 * 2 * 12 * 3600 * 1000L + 10000;

    public static PublicKey getPublicKey(String key) throws Exception {
        byte[] keyBytes;
        keyBytes = (new BASE64Decoder()).decodeBuffer(key);
        X509EncodedKeySpec keySpec = new X509EncodedKeySpec(keyBytes);
        KeyFactory keyFactory = KeyFactory.getInstance("RSA");
        PublicKey publicKey = keyFactory.generatePublic(keySpec);
        return publicKey;
    }

    /**
     * 公钥解析token
     *
     * @param token     用户请求中的token
     * @param publicKey 公钥
     * @return Jws<Claims>
     */
    private static Jws<Claims> parserToken(String token, PublicKey publicKey) {
        return Jwts.parser().setSigningKey(publicKey).parseClaimsJws(token);
    }

    public static boolean verify(String token)  {
        try {
            Resource resource = new ClassPathResource("public_key.txt");
            String publicKey = new String(FileCopyUtils.copyToByteArray(resource.getInputStream()));
            PublicKey pk = getPublicKey(publicKey);
            parserToken(token, pk);
        }catch (Exception ex){
            log.error(ex.getMessage());
            return false;
        }
        return true;
    }

    public static Jws<Claims> getClaims(String token){
        try {
            Resource resource = new ClassPathResource("public_key.txt");
            String publicKey = new String(FileCopyUtils.copyToByteArray(resource.getInputStream()));
            PublicKey pk = getPublicKey(publicKey);
            Jws<Claims> clasms = parserToken(token, pk);
            return clasms;
        }catch (Exception ex){
            log.error(ex.getMessage());
            return null;
        }
    }

    /**
     * 根据token获取用户id
     *
     * @param token 用户token
     * @return 用户id
     */
    public static long getUserIdByToken(String token){
        return new Long(getClaims(token).getBody().get("userId").toString());
    }

    /**
     * 获取登录用户id
     *
     * @return 用户id
     */
    public static long getUserId(){
        String token = (String) SecurityUtils.getSubject().getPrincipal();

        token = Optional.ofNullable(token).orElse("");
        if (StrUtil.isBlankOrUndefined(token)) {
            log.error("未获取到用户token");
            return 0L;
        }
        return JwtUtils.getUserIdByToken(token.toString());
    }

    public static long getTokenExp(String token){
        return new Long(getClaims(token).getBody().get("expire").toString());
    }

    public static JSONObject getPayload(String token){
        String payload = token.split("\\.")[1];
        String json = Base64.decodeStr(payload);
        return JSONUtil.parseObj(json);
    }
}
