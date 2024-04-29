package com.mgi.pacs.common.core.config;

import cn.hutool.core.collection.CollectionUtil;
import com.fasterxml.jackson.core.JsonGenerator;
import com.fasterxml.jackson.databind.JsonSerializer;
import com.fasterxml.jackson.databind.SerializerProvider;
import lombok.extern.slf4j.Slf4j;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

/**
 * 自定义Long集合序列化成String集合
 *
 * @author linzhunsheng
 * @date 2021/10/18
 */
@Slf4j
public class LongListJsonSerializer extends JsonSerializer<List<Long>> {

    @Override
    public void serialize(List<Long> value, JsonGenerator jsonGenerator, SerializerProvider serializerProvider) throws IOException {
        String text = "";
        //是否为空
        if (!CollectionUtil.isEmpty(value)) {
            try {
                List<String> strList = new ArrayList<>(value.size());
                for (Long aLong : value) {
                    strList.add(String.valueOf(aLong));
                }
                //格式化是否为空
                if (!CollectionUtil.isEmpty(strList)) {
                    jsonGenerator.writeObject(strList);
                    return;
                }
            } catch (Exception e) {
                log.error("解析长整形集合错误", e);
            }
        }
        jsonGenerator.writeString(text);
    }
}
