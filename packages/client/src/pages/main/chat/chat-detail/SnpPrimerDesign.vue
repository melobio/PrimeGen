<template>
    <div class="snp-primer-design" v-if="operations.length > 0">
        <template v-if="optionInfo?.stage">
            <v-form ref="formRef">
                <div class="option-item" v-for="(item, index) in operations" :key="item.title">
                    <span class="option-item-title">
                        <span class="required-icon" v-if="item?.type.includes('required')">*</span>
                        {{
                            item.key.split('_').join(' ').toUpperCase() }}
                    </span>
                    <div class="option-item-opt">
                        <v-text-field type="number" v-if="item.type.includes('input')" clearable
                            class="option-item-opt-input" v-model="operations[index].value[0]" :placeholder="`please input`"
                            :rules="getRules(item.key)" :disabled="isSubmited" hide-details="auto"
                            :label="item.title"></v-text-field>
                    </div>
                </div>
                <div style="display: flex;justify-content: flex-end;">
                    <v-btn variant="elevated" class='ml-3 mt-3' color='#2db489' @click="submitOption"
                        :disabled="!canSubmit">
                        {{ isSubmited ? `Submited` : `Submit` }}
                    </v-btn>
                </div>
            </v-form>
        </template>
    </div>
</template>
  
<script setup lang='ts'>
import { computed, ref, onMounted } from "vue";
import { useChatStore } from '@/stores/chats';
import { PrimerDesignFunctionType, type Message, type Operation, type SnpPrimerDesignInfo } from "@/stores/chat-types";
const chatStore = useChatStore();
const props = defineProps<Message>()
const operations = ref<Operation[]>([]);
const isSubmited = ref(false);
const formRef = ref()
// 1,Maximum amplicon length(最大放大长度)大于200小于500，
// 2,Minimum amplicon length(最小放大长度)大于100小于300
// 3,temperature(退火温度)30-80
// 4,Minimum GC content(最小GC含量)20-70
// 5,Maximum GC content(最大GC含量)20-70
// 6,Maximum length of a primer(最大引物长度)18-25
const optionRules = [
    {
        key: "min_amp_len",
        rules: [
            (value: any) => {
                if (value == '' || value <= 100 || value >= 300) {
                    return 'Minimum amplicon length greater than 100 less than 300'
                }
                if(getOperationValue("min_amp_len")>getOperationValue("max_amp_len")){
                    return 'Minimum amplicon length less than Maximum amplicon length'
                }
                return true

            }]
    },
    {
        key: "max_amp_len",
        rules: [
            (value: any) => {
                if (value == '' || value <= 200 || value >= 500) {
                    return 'Maximum amplicon length greater than 200 less than 500'
                }
                if(getOperationValue("max_amp_len")<getOperationValue("min_amp_len")){
                    return 'Maximum amplicon length greater than Minimum amplicon length'
                }
                return true;
            }]
    },
    {
        key: "temperature",
        rules: [
            (value: any) => {
                if (value == '' || value <= 30 || value >= 80) {
                    return 'temperature greater than 30 less than 80'
                }
                return true
            }]
    },
    {
        key: "min_GC",
        rules: [
            (value: any) => {
                if (value == '' || value <= 0 || value >= 1) {
                    return 'Minimum GC content greater than 0 less than 1'
                }
                if(getOperationValue("min_GC")>=getOperationValue("max_GC")){
                    return 'Minimum GC content less than Maximum GC content'
                }
                return true
            }]
    },
    {
        key: "max_GC",
        rules: [
            (value: any) => {
                if (value == '' || value <= 0 || value >= 1) {
                    return 'Maximum length of a primer greater than 0 less than 1'
                }
                if(getOperationValue("max_GC")<getOperationValue("min_GC")){
                    return 'Maximum GC content greater than Minimum GC content'
                }
                return true
            }]
    },
    {
        key: "max_primer_len",
        rules: [
            (value: any) => {
                if (value == '' || value <= 18 || value >= 25) {
                    return 'Maximum GC content greater than 18 less than 25'
                }
                return true
            }]
    }
]

onMounted(() => {
    operations.value = props?.optionInfo?.operations || [];
    isSubmited.value = props?.optionInfo?.submitted ? true : false
})

const validate = async () => {
    const { valid } = await formRef.value.validate()
    return valid
}

const getOperationValue = (key:string)=>{
    let targetValue = '';
    operations.value.forEach(operation=>{
        if(operation.key == key){
            targetValue = operation.value[0];
        }
    })
    return targetValue;
}
const getRules = (key: string) => {
    return optionRules.find(rule => rule?.key == key)?.rules || []
}

const canSubmit = computed(() => {
    return (!chatStore.currentChat.isGenerating) && !isSubmited.value;
})

const submitOption = async () => {
    if (await validate()) {
        isSubmited.value = true;
        let stage = props.optionInfo?.stage ?? 2;
        let userReply = "I have submitted";
        let state = props.optionInfo?.state || 'continue';
        let primer_design_prompt = props.optionInfo?.primer_design_prompt || '';
        let primer_design_dict = props.optionInfo?.primer_design_dict || {};
        let primer_type = props.optionInfo?.primer_type || PrimerDesignFunctionType.snp_primer_design_type;
        operations.value.forEach((item, index) => {
            if (!Array.isArray(item.value)) {
                operations.value[index] = { ...item, value: [item.value] }
            }
        })
        let data: SnpPrimerDesignInfo = {
            stage,
            state,
            primer_type,
            primer_design_prompt,
            primer_design_dict,
            operations: operations.value,
            // 上一次返回的值
            data: props?.optionInfo?.data ?? {},
        }
        await chatStore.sendSelectOptionToPrimerDesign(userReply, data, props)
    }
}
</script>
  
<style scoped lang='scss'>
.snp-primer-design {
    background-color: #1E2032;
    border-radius: 10px;
    padding: 10px 20px 10px;

    .option-item {
        width: 100%;
        margin-top: 20px;

        &-title {
            font-size: 14px;
            font-weight: 600;
            position: relative;
        }

        &-des {
            position: relative;
        }

        &-opt {
            margin-top: 10px;

            &-input {}

        }

        .required-icon {
            color: #dc362e;
            margin-right: 2px;
            position: absolute;
            left: -5px;
            top: 0px;
        }
    }

    .title {
        color: #ECECEC;
        font-size: 18px;
        line-height: 44px;
        height: 44px;
    }

}
</style>