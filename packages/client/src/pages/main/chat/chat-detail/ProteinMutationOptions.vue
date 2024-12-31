<template>
    <div class="protein-mutation-options" v-if="operations.length">
        <!-- let user choose -->
        <template v-if="props?.optionInfo?.stage == 2">
            <v-form ref="formRef">
                <div class="option-item" v-for="(item, index) in operations" :key="item.title">
                    <span class="option-item-title">
                        <span class="required-icon" v-if="item?.type.includes('required')">*</span>
                        {{
                            item.key.split('_').join(' ').toUpperCase() }}
                    </span>
                    <div class="option-item-opt">
                        <v-text-field type="string" v-if="item.type.includes('input')" clearable
                            class="option-item-opt-input" v-model="operations[index].value[0]" :placeholder="`please input`"
                            :disabled="isSubmited" hide-details="auto"
                            :label="item.title"></v-text-field>
                    </div>
                </div>
                <div style="display: flex;justify-content: flex-end;">
                    <v-btn variant="elevated" class='ml-3 mt-3' color='#2db489' @click="submitOption"
                        :disabled="chatStore.currentChat.isGenerating">
                        {{ isSubmited ? `Submited` : `Submit` }}
                    </v-btn>
                </div>
            </v-form>
        </template>
        <template v-else-if="props?.optionInfo?.stage == 5">
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
                            :rules="getRules(item.key)"
                            :disabled="isSubmited" hide-details="auto"
                            :label="item.title"></v-text-field>
                    </div>
                </div>
                <div style="display: flex;justify-content: flex-end;">
                    <v-btn variant="elevated" class='ml-3 mt-3' color='#2db489' @click="submitOption"
                        :disabled="!lengthCanSubmit">
                        {{ isSubmited ? `Submited` : `Submit` }}
                    </v-btn>
                </div>
            </v-form>
        </template>
        <template v-else >
            <div v-for="(operation, index) in operations" :key="index" class="option-item">
                <v-table style="background: transparent;"
                    v-if="operation.type.includes('single') || operation.type.includes('multi')">
                    <thead>
                        <tr>
                            <th class="text-center">
                                Protein Name
                            </th>
                            <th class="text-center">
                                Protein Length
                            </th>
                            <th class="text-center">
                                Organism
                            </th>
                            <th class="text-center">
                                Function Information
                            </th>
                            <th class="text-center">
                                Select
                            </th>
                        </tr>
                    </thead>
                    <tbody>
                        <template v-if="operations.length">
                            <template v-for="(option) in operations">
                                <tr class="text-center" v-for="(item, index) in option.options" :key="index">
                                    <td class="text-center">{{ item['Protein Name'] }}</td>
                                    <td class="text-center">{{ item['Length'] }}</td>
                                    <td class="text-center">{{ item['Organism'] }}</td>
                                    <td class="text-center">
                                        <v-btn density="compact" icon="" variant="plain">
                                            <span class="mdi mdi-dots-horizontal"></span>
                                            <v-tooltip activator="parent" location="end">
                                                <div style="width:400px;">
                                                    {{ item['Function information'] }}
                                                </div>
                                            </v-tooltip>
                                        </v-btn>
                                    </td>
                                    <td class="text-center">
                                        <v-checkbox hide-details :disabled="isSubmited"
                                            @change="() => selectedChange(index)" v-model="selected"
                                            :value="index"></v-checkbox>
                                    </td>
                                </tr>
                            </template>
                        </template>
                        <tr v-else>
                            <td class="text-center">
                                -
                            </td>
                            <td class="text-center">
                                -
                            </td>
                            <td class="text-center">
                                -
                            </td>
                            <td class="text-center">
                                -
                            </td>
                            <td class="text-center">
                                -
                            </td>
                        </tr>
                    </tbody>
                </v-table>
                <template v-else>
                    <span class="option-item-title">
                        {{ operation.title }}
                    </span>
                    <template v-if="uploadedFiles.length > 0 && isSubmited">
                        <div class="file-uploaded" v-for="(fileItem, fileIndex) in uploadedFiles" :key="fileIndex">
                            {{ fileItem.split('/').pop() }}
                        </div>
                    </template>
                    <template v-else>
                        <v-file-input @change="handleFileChange" v-model="selectedFiles" prepend-icon=""
                            label="Upload Gene File" :disabled="isSubmited" :rules="rules">
                        </v-file-input>
                    </template>
                </template>
                <div style="display: flex;justify-content: flex-end;">
                    <v-btn :disabled="!canSubmit" variant="elevated" class='ml-3 mt-3' color='#2db489' @click="submitOption">
                        {{ isSubmited ? `Submited` : `Submit` }}
                    </v-btn>
                </div>
            </div>
        </template>
    </div>
</template>
<script lang="ts" setup>
import { computed, ref, onMounted } from "vue";
import { useChatStore } from '@/stores/chats';
import { NCBIFunctionType, type Message, type ProteinMutationLength, type ProteinMutationInfo, type Operation } from "@/stores/chat-types";
import { AgentType } from "@xpcr/common/src/agent";
import { uploadSequenceFiles } from '@/api';

const selectedFiles = ref<File[]>([]);
const uploadedFiles = ref<string[]>([]);
const max_length = ref(0);

const fileExtensionRegex = /\.(fasta|fastq|fna|fa|csv|tsv)$/i;

const chatStore = useChatStore();
const props = defineProps<Message>();
const selected = ref<Array<any>>([]);
const isSubmited = ref(false);
const operations = ref<Operation[]>([]);
const rules = ref([
    (files: any) => {
        let validate = true;
        files.forEach((file: any) => {
            if (!fileExtensionRegex.test(file.name)) {
                validate = false;
            }
        });
        return validate || 'You can only upload fasta/fastq/fna/fa/csv/tsv files';
    },
],);


// 1. start position
// 2. end position
const optionRules = [
    {
        key: "start_pos",
        rules: [
            (value: any) => {
                if (value == '' || value < 1) {
                    return 'Start position less than 1'
                }
                if(Number(getOperationValue("start_pos")) >= Number(getOperationValue("end_pos"))){
                    return 'Start position larger than end position'
                }
                if (Number(getOperationValue("end_pos")) - Number(getOperationValue("start_pos")) > 30000) {
                    return 'End position - Start position should not larger than 30000!'
                }
                return true

            }]
    },
    {
        key: "end_pos",
        rules: [
            (value: any) => {
                if (value == '' || value > max_length.value) {
                    return 'End position larger than ' + max_length.value
                }
                if(Number(getOperationValue("start_pos")) >= Number(getOperationValue("end_pos"))){
                    return 'Start position larger than end position'
                }
                if (Number(getOperationValue("end_pos")) - Number(getOperationValue("start_pos")) > 30000) {
                    return 'End position - Start position should not larger than 30000!'
                }
                return true;
            }]
    }
]

const getOperationValue = (key:string)=>{
    let targetValue = 0;
    operations.value.forEach(operation=>{
        if(operation.key == key){
            targetValue = operation.value[0];
        }
    })
    return targetValue;
}

const getRules = (key: string) => {
    return optionRules.find(rule => rule?.key == key)?.rules || [];
}

onMounted(() => {
    operations.value = props.optionInfo?.operations || [];
    isSubmited.value = props.optionInfo?.submitted ?? false;
    max_length.value = props.optionInfo?.data?.max_length ?? 0;

    operations.value.forEach(operation => {
        if (!operation.type.includes('file')) {
            operation.options?.forEach((item: Record<string, any>, index: number) => {
                operation.value.forEach(operationItem => {
                    if (operationItem['Protein ID'] == item['Protein ID']) {
                        selected.value.push(index)
                    }
                })
            })
        } else {
            if (operation.value.length > 0) {
                uploadedFiles.value = [...operation.value];
            }
        }
    });
})

const canSubmit = computed(() => {
    return (!chatStore.currentChat.isGenerating) && (selected.value?.length > 0 || selectedFiles.value.length>0)&& !(isSubmited.value);
})

const lengthCanSubmit = computed(() => {
    return (!chatStore.currentChat.isGenerating) && !isSubmited.value && operations.value[0].value[0] > 0 && operations.value[1].value[0] > 0;
})

// single selection
const selectedChange = (index: number) => {
    selected.value = [index]
}

const submitOption = async () => {
    isSubmited.value = true;
    let stage = props.optionInfo?.stage ?? 2;
    let protein_mutation_dict = props.optionInfo?.protein_mutation_dict || {};
    let search_type = props.optionInfo?.search_type ?? NCBIFunctionType.protein_mutation_type;
    let userReply = "I have submitted";

    operations.value.forEach((item, index) => {
        if (!item.type.includes('file') && selected.value[0]) {
            operations.value[index] = { ...item, value: [item.options[selected.value[0]]] }
        }
    })
    let data: ProteinMutationInfo = { stage, operations: operations.value, protein_mutation_dict, search_type, data: props?.optionInfo?.data ?? {} }
    await chatStore.sendSelectOptionToNcbiSearch(userReply, data, props)

}
const handleFileChange = () => {
    let formData = new FormData();
    formData.append('files', selectedFiles.value[0]);
    uploadSequenceFiles(formData).then(({ data }) => {
        if (data?.success == 'True') {
            operations.value.forEach((item, index) => {
                if (item.type.includes('file')) {
                    operations.value[index].value = data.filenames
                }
            })
        }
    });
}
</script>
<style scoped lang="scss">
.protein-mutation-options {
    padding: 10px 20px;
    box-sizing: border-box;
    border-radius: 24px;
    border: 3px solid transparent;
    background: #151728;
    margin-top: 20px;

    .option-item {
        width: 100%;
        margin-top: 20px;

        &:first-child {
            margin-top: 0px;
        }

        &-title {
            font-size: 14px;
            font-weight: 600;
            margin-bottom: 10px;
            display: block;
        }

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
</style>
