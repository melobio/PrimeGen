package com.mgi.pacs.primer.constant;

public enum CureKeysFromEnum {

    //定义枚举对象 static final
    FROM_BGI_FAMILY("bgi_family", 1);

    private final String name;
    private final Integer type;

    CureKeysFromEnum(String name, Integer stateNum) {
        this.name = name;
        this.type = stateNum;
    }

    public Integer getType() {
        return type;
    }

    public String getName() {
        return name;
    }

    //toString默认是变量名，

    public static CureKeysFromEnum getTypeByName(String name) {
        for (CureKeysFromEnum item : CureKeysFromEnum.values()) {
            if (item.name.equals(name)) {
                return item;
            }
        }
        throw new IllegalArgumentException("CureKeysFrom状态参数异常");
    }
}
