package com.mgi.pacs.common.core.config;

import com.fasterxml.jackson.core.JsonParser;
import com.fasterxml.jackson.core.type.TypeReference;
import com.fasterxml.jackson.databind.DeserializationContext;
import com.fasterxml.jackson.databind.JsonDeserializer;
import com.fasterxml.jackson.databind.ObjectMapper;
import lombok.extern.slf4j.Slf4j;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

/**
 * 自定义String集合反序列化成Long集合
 *
 * @author linzhunsheng
 * @date 2021/11/11
 */
@Slf4j
public class LongListJsonDeserializer extends JsonDeserializer<List<Long>> {

    static final TypeReference<List<Long>> LIST_TYPE = new TypeReference<List<Long>>() {};

    @Override
    public List<Long> deserialize(JsonParser jsonParser, DeserializationContext deserializationContext) throws IOException {
        ObjectMapper mapper = new ObjectMapper();
        List<Long> longList = new ArrayList<>();
        try {
            longList = mapper.readValue(jsonParser, LIST_TYPE);
        } catch (IOException e) {
            log.error("反序列化长整形集合错误", e);
        }
        return longList;
    }
}
