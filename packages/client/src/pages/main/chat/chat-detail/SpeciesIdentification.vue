<template>
    <template v-if="optionInfo?.stage == 1">
        <v-form ref="formRef">
            <div class="species-identification-options" >
                <div v-for="(operation, index) in operations" :key="index" class="option-item">
                    <v-table v-if="operation?.widget == 'table' && (operation.type.includes('single') || operation.type.includes('multi'))" style="background: transparent;">
                        <thead>
                        <tr>
                            <th class="text-center">Select</th>
                            <th v-for="head in getTableHeads()" :key="head" class="text-center">
                            {{ formatColumnName(head) }}
                            </th>
                        </tr>
                        </thead>
                        <tbody>
                            <tr v-for="(item, itemIndex) in operation.options" :key="itemIndex" class="text-center">
                                <td class="text-center">
                                    <v-checkbox
                                        hide-details
                                        :disabled="isSubmited"
                                        @change="() => selectedChange(itemIndex, operation.type)"
                                        v-model="selected"
                                        :value="itemIndex"
                                    ></v-checkbox>
                                </td>
                                <td v-for="head in getTableHeads()" :key="head" class="text-center">
                                    <template v-if="head.toLowerCase() == 'ftp path'">
                                        <a :href="item[head.replace(' ', '_').toLowerCase()]" download target="_blank" rel="noopener noreferrer">Detail</a>
                                    </template>
                                    <template v-else >
                                        {{ item[head.replace(' ', '_').toLowerCase()] }}
                                    </template>
                                </td>
                            </tr>
                        </tbody>
                    </v-table>
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
    <template v-else-if="operations.length > 0">
        <div class="species-identification" >
            <template v-if="optionInfo?.stage">
                <v-form ref="formRef">
                    <div class="option-item" v-for="(item, index) in operations" :key="item.title">
                        <span class="option-item-title">
                            <span class="required-icon" v-if="item?.type.includes('required')">*</span> {{
                                item.key.split('_').join(' ').toUpperCase() }}
                        </span>
                        <div class="option-item-des">
                            {{ item.title }} :
                        </div>
                        <div class="option-item-opt">
                            <!-- 渲染选择题(单选或多选) -->
                            <div v-if="item.type.includes('single') || item.type.includes('multi')"
                                :key="JSON.stringify(ifSelectedMicrobe)">
                                <template v-if="isSubmited">
                                    <div class="file-uploaded">
                                        <div v-for="(text, index) in item.value" :key="text">
                                            {{ text }}<span v-if="item.value.length - 1 > index">,</span>
                                        </div>
                                    </div>
                                </template>
                                <template v-else>
                                    <v-combobox v-if="item.type.includes('input')" clearable class="mx-auto"
                                        :multiple="ifMultiple(item)" v-model="operations[index].value" 
                                        :rules="item?.type.includes('required')?optionRules:[]"
                                        :disabled="isSubmited" :placeholder="`please select`"
                                        :items="item.options"></v-combobox>
                                    <v-select v-else clearable class="option-item-opt-select" v-model="operations[index].value"
                                        :multiple="ifMultiple(item)" :placeholder="`please select`" 
                                        :rules="item?.type.includes('required')?optionRules:[]"
                                        :disabled="isSubmited" :items="item.options"></v-select>
                                </template>
                            </div>
                            <!-- 渲染文件上传 -->
                            <template v-if="item.type.includes('file')">
                                <template v-if="item.value.length > 0 && isSubmited">
                                    <div class="file-uploaded" v-for="(fileItem, fileIndex) in operations[index].value"
                                        :key="fileIndex">
                                        {{ fileItem.split('/').pop() }}
                                    </div>
                                </template>
                                <template v-else>
                                    <v-file-input class="file-input" v-model="selectedFiles[item.key]"
                                        @change="handleFileChange(item.key)" :label="`please upload`" prepend-icon=""
                                        :disabled="isSubmited" :rules="fileRules">
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
</template>
<script lang="ts" setup>
import { computed, ref, onMounted, watch } from "vue";
import { useChatStore } from '@/stores/chats';
import { NCBIFunctionType, type Message, type SpeciesIdentificationInfo, type Operation } from "@/stores/chat-types";
import { uploadSequenceFiles } from "@/api";
const chatStore = useChatStore();
const props = defineProps<Message>();
const selectedFiles = ref<Record<string, any>>({});
// 打勾选项
const selected = ref<Array<any>>([]);
const operations = ref<Operation[]>([]);
const formRef = ref()
const optionRules = [
    (value: any) => {
        if (value&&value.length!==0) {return true}
        return 'must be selected.'
    },
]

// 定义正则表达式匹配菌株名_cds_from_genomic.fna的模式
// *_gene.fna  // target-gene non-gene

const strainFileRegex = /^(.+?)_cds_from_genomic\.(fna|fa)$/;
const geneFileregex = /^(.+?).(fna|fa|fasta)$/;
const validateRules = [
    {
        keys: ['strain_path'],
        regexRule: strainFileRegex,
        tips: 'The file you uploaded is incorrect. We only support files containing coding sequences, i.e., CDS files. Additionally, please rename your file as follows: "strain_name_cds_from_genomic.fna".'
    },
    {
        keys: ['non_target_gene_path', 'target_gene_path'],
        regexRule: geneFileregex,
        tips: 'The file you upload does not meet the prescribed format.Please make sure the file suffix is ”*.fa“ or "*.fna" or "*.fasta"'

    }
]
const validate = async () => {
    const { valid } = await formRef.value.validate()
    return valid
}

const fileRules = ref([
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

watch(() => operations.value[0]?.value, (newValue, oldValue) => {
    if (!isSubmited.value && operations.value[0]?.key == 'experiment_select') {
        operations.value.forEach((item, index) => {
            if (item.key == 'strain_select') {
                operations.value[index] = { ...item, value: [] };
            }
        })
    }
}, {
    deep: true,
})
const isSubmited = ref(false);

const ifSelectedMicrobe = computed(() => {
    return operations.value.filter(item => item.key == 'experiment_select' && item.value?.includes('microbe')).length > 0

})

const ifMultiple = (operation: Operation) => {
    return ifSelectedMicrobe.value ? false : operation.type.includes('multi')
}

const getTableHeads = () => {
    return ['organism name', 'group', 'assembly level', 'representative', 'genus', 
            'species', 'infraspecific name', 'isolate', 'accession', 'taxid', 'origin', 'ftp path']
}

// 转换成首字母大写
const formatColumnName = (name: string) => {
  return name.replace(/(?:^|\s)\w/g, match => match.toUpperCase());
}

// 勾选了就往勾选列表增加offset
const selectedChange = (index: number, operation_type: Array<string>) => {
    if (operation_type.includes('single')) {
        selected.value = [index]
    } else {
        selected.value.push(index)
    }
}

onMounted(() => {
    operations.value = props?.optionInfo?.operations || [];
    isSubmited.value = props?.optionInfo?.submitted ? true : false
    operations.value.forEach(item => {
        if (item && item?.type.includes('file')) {
            selectedFiles.value[item.key] = '';
        }
    })
})

const canSubmit = computed(() => {
    return (!chatStore.currentChat.isGenerating) && !isSubmited.value;
})

watch(() => selectedFiles.value, (newVal, oldVal) => {
    for (const key in selectedFiles.value) {
        const files = selectedFiles.value[key];
        operations.value.forEach((item, index) => {
            if (item.type.includes('file') && item.key == key) {
                operations.value[index].value = files
            }
        })
    }
}, { deep: true })


const handleFileChange = (key: string) => {
    let formData = new FormData();
    formData.append('files', selectedFiles.value[key][0]);
    uploadSequenceFiles(formData).then(({ data }) => {
        if (data?.success == 'True') {
            operations.value.forEach((item, index) => {
                if (item.type.includes('file') && item.key == key) {
                    operations.value[index].value = data.filenames
                }
            })
        }
    });
}
const submitOption = async () => {
    if (await validate()) {
        isSubmited.value = true;
        let stage = props.optionInfo?.stage ?? 2;
        let userReply = "I have submitted";
        let upload_file_flag = false;

        let species_identification_dict = props.optionInfo?.species_identification_dict || {}
        let search_type = props.optionInfo?.search_type || NCBIFunctionType.species_identification_type;
        operations.value.forEach((item, index) => {
            if (!Array.isArray(item.value)) {
                operations.value[index] = { ...item, value: [item.value] }
            }
            if (item.type.includes('file')&&item.value.length>0) {
                upload_file_flag = true;
            }
            // 非文件类的, 选项
            else if (selected.value) {
                operations.value[index] = { ...item, value: [item.options[selected.value[0]]] }
            }
        })
        let data: SpeciesIdentificationInfo = {
            stage,
            search_type,
            species_identification_dict,
            operations: operations.value,
            upload_file_flag,
            // 上一次返回的值
            data: props?.optionInfo?.data ?? {},
        }
        await chatStore.sendSelectOptionToNcbiSearch(userReply, data, props)
    }
}
</script>
<style scoped lang="scss">
.species-identification-options {
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
            position: relative;
        }

        &-des {
            position: relative;
        }

        &-opt {
            margin-top: 10px;

            &-select {}

        }

        .required-icon {
            color: #dc362e;
            margin-right: 2px;
            position: absolute;
            left: -5px;
            top: 0px;
        }
    }

    .file-input {
        .v-messages {
            .v-messages__message {
                padding: 5px 0;
            }
        }
    }

    .file-uploaded {
        border-radius: 3px;
        min-height: 50px;
        padding: 0 12px;
        border: 1px solid #ccc;
        color: #ccc;
        pointer-events: none;
        opacity: 0.4;
        display: flex;
        flex-direction: column;
        justify-content: center;
    }
}</style>
