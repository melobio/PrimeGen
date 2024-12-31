<template>
    <div v-if="operations.length > 0" class="protocol-select-info">
        <v-form ref="formRef">
            <!-- user choose Protocol Library -->
            <template v-if="Number(stage) == 2" >
                <template v-for="operation in operations" :key="operation.key" >

                    <template v-if="(operation.type.includes('select'))">
                        <div class="option-item">
                            <span class="option-item-title">
                                {{ operation.title }}
                            </span>
                            <span>
                                <v-select size="large" class="option-item-opt-select" v-model="operation.value[0]" clearable
                                    :rules="optionRules" :disabled="isSubmited" :items="protocolSelections">
                                </v-select>
                            </span>
                        </div>
                    </template>

                    <template v-if="(operation.type.includes('file'))">
                        <div class="option-item">
                            <span class="option-item-title">
                                {{ operation.title }}
                            </span>
                            <!-- uploaded -->
                            <template v-if="uploadedFiles.length > 0 && isSubmited">
                                <div class="file-uploaded" v-for="(fileItem, fileIndex) in uploadedFiles" :key="fileIndex">
                                    {{ fileItem.split('/').pop() }}
                                </div>
                            </template>
                            <!-- wait to upload -->
                            <template v-else>
                                <v-file-input @change="handleFileChange" v-model="selectedFiles" prepend-icon=""
                                    label="Upload Protocol File" :disabled="isSubmited" :rules="fileRules">
                                </v-file-input>
                            </template>
                        </div>
                    </template>

                </template>
            </template>

            <!-- step 2: according to the Protocol Library user selected, input parameters -->
            <!-- when protocol_design_type == 'auto_protocol_design' -->
            <template v-if="Number(stage) == 3 || Number(stage) == 4" >
                <template v-for="param in operations" :key="param.key" >
                    <div class="option-item">
                        <span class="option-item-title">
                            {{ param.title }}
                        </span>
                        <div class="option-item-opt">
                            <v-text-field
                                class="option-item-opt-input"
                                type="string"
                                size="large"
                                v-model="param.value[0]"
                                :label="param.title"
                                :disabled="isSubmited"
                                :rules="optionRules">
                            </v-text-field>
                        </div>
                    </div>
                </template>
            </template>
            <template v-if="props?.optionInfo?.state !== 'stop'">
                <div style="display: flex;justify-content: flex-end;">
                    <v-btn :disabled="!canSubmit" variant="elevated" class='ml-3 mt-3' color='#2db489' @click="submitOption">
                        {{ isSubmited ? `Submitted` : `Submit` }}
                    </v-btn>
                </div>
            </template>

        </v-form>
    </div>
</template>

<script lang="ts" setup>
import { computed, ref, onMounted } from "vue";
import { useChatStore } from '@/stores/chats';
import { type Message, type Operation, type ProtocolDesignInfo, type ProtocolDesignResultInfo, ProtocolDesignFunctionType } from "@/stores/chat-types";
import { uploadSequenceFiles } from '@/api';


const chatStore = useChatStore();

const props = defineProps<Message>();

const optionInfo = ref<ProtocolDesignResultInfo>();

const operations = ref<Operation[]>([]);

const stage = ref(1);

const isSubmited = ref(false);

const uploadedFiles = ref<string[]>([]);

const selectedFiles = ref<File[]>([]);

const fileExtensionRegex = /\.(pdf|txt|doc|docx)$/i;

const fileRules = ref([
    (files: any) => {
        let validate = true;
        let tips = '';
        files.forEach((file: any) => {
            if (!fileExtensionRegex.test(file.name)) {
                validate = false;
                tips = 'File type not support, only PDF, TXT, DOC and DOCX are supported';
            }
        });
        return validate ? true : tips
    }
]);

const protocolSelections = ref<Array<string>>([]);

const optionRules = [
    (value: any) => {
        if (value) return true
        return 'must be select or input parameter'
    },
]


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

onMounted(() => {
    optionInfo.value = props?.optionInfo;
    operations.value = props?.optionInfo?.operations || [];
    isSubmited.value = props?.optionInfo?.submitted ? true:false
    stage.value = props?.optionInfo?.stage || 1;
    operations.value.forEach(operation => {
        if (operation.type.includes('file')) {
            uploadedFiles.value = [...operation.value];
        }
        else {
            protocolSelections.value = [...operation.options];
        }
    });
})

const canSubmit = computed(() => {
    const isModified = operations.value.some(operation => operation.value.length > 0);
    return !chatStore.currentChat.isGenerating && !isSubmited.value && isModified;
})

const submitOption = async () => {
    isSubmited.value = true;
    let userReply = "I have submitted";

    let data: ProtocolDesignInfo = {
        stage: props.optionInfo?.stage ?? 2,
        protocol_design_type: props?.optionInfo?.protocol_design_type ?? ProtocolDesignFunctionType.template_protocol_design,
        operations: operations.value,
        data: props?.optionInfo?.data ?? {},
    }
    await chatStore.sendOptionToProtocolDesign(userReply, data, props)

}


</script>

<style scoped lang="scss">
.protocol-select-info {
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
