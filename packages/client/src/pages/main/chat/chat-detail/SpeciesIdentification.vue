<template>
    <div>
        <template v-if="optionInfo?.stage">
            <v-form v-if="operations.length > 0" class="species-identification">
                <div class="option-item" v-for="(item, index) in operations" :key="item.title">
                    <span class="option-item-title">
                        {{ item.title.toUpperCase() }}:
                    </span>
                    <div class="option-item-opt">
                        <!-- 渲染选择题(单选或多选) -->
                        <template v-if="item.type.includes('single') || item.type.includes('multi')">
                            <v-combobox v-if="item.type.includes('input')" clearable class="mx-auto"
                                :multiple="item.type.includes('multi')" v-model="operations[index].value"
                                :rules="optionRules" :disabled="isSubmited" :placeholder="`Select ${item.title}`"
                                :items="item.options"></v-combobox>
                            <v-select v-else clearable class="option-item-opt-select" v-model="operations[index].value"
                                :multiple="item.type.includes('multi')" :placeholder="`Select ${item.title}`"
                                :rules="optionRules" :disabled="isSubmited" :items="item.options"></v-select>
                        </template>
                        <!-- 渲染文件上传 -->
                        <template v-if="item.type.includes('file')">
                            <template v-if="item.value.length > 0 && isSubmited">
                                <div class="file-uploaded" v-for="(fileItem, fileIndex) in operations[index].value"
                                    :key="fileIndex">
                                    {{ fileItem.split('/').pop() }}
                                </div>
                            </template>
                            <template v-else>
                                <v-file-input v-model="selectedFiles[item.key]" @change="handleFileChange(item.key)"
                                    :label="item.title" prepend-icon="" :disabled="isSubmited" :rules="rules">
                                </v-file-input>
                            </template>
                        </template>
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
<script lang="ts" setup>
import { computed, ref, onMounted } from "vue";
import { useChatStore } from '@/stores/chats';
import { NCBIFunctionType, type Message, type SpeciesIdentificationInfo, type Operation } from "@/stores/chat-types";
import { AgentType } from "@xpcr/common/src/agent";
import { uploadSequenceFiles } from "@/api";
const chatStore = useChatStore();
const props = defineProps<Message>();
const selectedFiles = ref<Record<string, any>>({});
const operations = ref<Operation[]>([]);
const optionRules = [
    (value: any) => {
        if (value) return true
        return 'must be selected.'
    },
]

// 定义正则表达式匹配菌株名_cds_from_genomic.fna的模式
// *_gene.fna  // target-gene non-gene

const strainFileRegex = /^(.+?)_cds_from_genomic\.(fna|fa)$/;
const geneFileregex = /^(.+?)_gene\.(fna|fa)$/;
const validateRules = [
    {
        keys: ['strain_path'],
        regexRule: strainFileRegex,
        tips: 'The file you upload does not meet the prescribed format.Please make sure the file name follows the format of "*_cds_from_genomic" and file suffix is fa or fna'
    },
    {
        keys: ['non_target_gene_path', 'target_gene_path'],
        regexRule: geneFileregex,
        tips: 'The file you upload does not meet the prescribed format.Please make sure the file name follows the format of "*_gene" and file suffix is fa or fna'

    }
]
const rules = ref([
    (files: any) => {
        let validateName = true;
        let validateTips = '';
        for (const key in files) {
            const file = files[key];
            operations.value.forEach(item => {
                if (item.type.includes('file')) {
                    let fileValidate = validateRules.find(rule => rule.keys.includes(item.key));
                    if (fileValidate && !fileValidate?.regexRule.test(file.name)) {
                        validateName = false;
                        validateTips = fileValidate.tips
                    }
                }
            })
        }
        let validateNameText = validateName || validateTips;
        return validateNameText;
    },
],);

const isSubmited = ref(false);
const upload_file_flag = ref(false);

const finishedUploadFile = computed(() => {
    let validate = true;
    for (const key in selectedFiles.value) {
        const files = selectedFiles.value[key];
        if (files && files.length > 0) {
            let fileValidate = validateRules.find(rule => rule.keys.includes(key));
            files.forEach((file: File) => {
                if (fileValidate && !fileValidate?.regexRule.test(file.name)) {
                    validate = false;
                }
            })
        }

    }

    return validate
})
onMounted(() => {
    operations.value = props?.optionInfo?.operations || [];
    operations.value.forEach(item => {
        if (item.value.length > 0) {
            isSubmited.value = true;
        }
        if (item.type.includes('file')) {
            selectedFiles.value[item.key] = '';
        }
    })
})

const canSubmit = computed(() => {
    // const lastSearchMsgId = chatStore.currentChatMessages.filter(item => item.agentType == AgentType.SEQUENCE_SEARCH).pop()?.id || '';
    // props.id == lastSearchMsgId && 
    let validateFormData = false;
    operations.value.forEach(operation=>{
        if(operation.key == 'experiment_select'){
            validateFormData = operation.value.length>0;
        }
    })
    
    return (!chatStore.currentChat.isGenerating) && !isSubmited.value && finishedUploadFile.value && validateFormData;
})

const handleFileChange = (key: string) => {
    let formData = new FormData();
    formData.append('files', selectedFiles.value[key][0]);
    uploadSequenceFiles(formData).then(({ data }) => {
        if (data?.success == 'True') {
            upload_file_flag.value = true;
            operations.value.forEach((item, index) => {
                if (item.type.includes('file') && item.key == key) {
                    operations.value[index].value = data.filenames
                }
            })
        }
    });
}
const submitOption = async () => {
    isSubmited.value = true;
    let stage = props.optionInfo?.stage ?? 2;
    let userReply = "I have submitted";
    let species_identification_dict = props.optionInfo?.species_identification_dict || {}
    let search_type = props.optionInfo?.search_type || NCBIFunctionType.species_identification_type;
    operations.value.forEach((item, index) => {
        if (!Array.isArray(item.value)) {
            operations.value[index] = { ...item, value: [item.value] }
        }
    })
    let data: SpeciesIdentificationInfo = {
        stage,
        search_type,
        species_identification_dict,
        operations: operations.value,
        upload_file_flag: upload_file_flag.value
    }
    await chatStore.sendSelectOptionToNcbiSearch(userReply, data, props)
}
</script>
<style scoped lang="scss">
.species-identification {
    padding: 10px 20px;
    box-sizing: border-box;
    border-radius: 24px;
    border: 3px solid transparent;
    background: #151728;
    margin: 10px 0;
    margin-top: 20px;

    .option-item {
        width: 100%;
        margin-top: 20px;

        &-title {
            font-size: 14px;
            font-weight: 600;
        }

        &-opt {
            margin-top: 10px;

            &-select {}

            .file-uploaded {
                border-radius: 3px;
                height: 50px;
                padding: 0 12px;
                line-height: 50px;
                border: 1px solid #ccc;
                color: #ccc;
                pointer-events: none;
                opacity: 0.4;
            }
        }
    }
}
</style>
