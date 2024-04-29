package com.mgi.pacs.primer.mapper;

import org.apache.commons.collections4.map.LinkedMap;
import org.apache.ibatis.annotations.Select;
import org.springframework.stereotype.Repository;

import java.util.List;
import java.util.Map;

@Repository
public interface CustomMapper {

    @Select("${sql}")
    List<LinkedMap<String, Object>> executeSelect(String sql);
}
