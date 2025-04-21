<template>

    <v-form ref="formRef">
        <!-- 步骤一: 用户选择Protocol Library, 选择完后上传文件 -->
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
                        <!-- 已经上传文件 / 已经提交 -->
                        <template v-if="uploadedFiles.length > 0 && isSubmited">
                            <div class="file-uploaded" v-for="(fileItem, fileIndex) in uploadedFiles" :key="fileIndex">
                                {{ fileItem.split('/').pop() }}
                            </div>
                        </template>
                        <!-- 等待上传文件 -->
                        <template v-else>
                            <v-file-input @change="handleFileChange" v-model="selectedFiles" prepend-icon=""
                                label="Upload Protocol File" :disabled="isSubmited" :rules="fileRules">
                            </v-file-input>
                        </template>
                    </div>
                </template>

            </template>
        </template>

        <!-- 步骤二: 根据用户选择的Protocol Library, 输入函数执行的参数信息 -->
        <!-- 步骤三: 当protocol_design_type == 'auto_protocol_design'时, 显示每个步骤需要输入的参数 -->
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
        <!-- 提交按钮 -->
        <template v-if="props?.optionInfo?.state !== 'stop'">
            <div style="display: flex;justify-content: flex-end;">
                <v-btn :disabled="!canSubmit" variant="elevated" class='ml-3 mt-3' color='#2db489' @click="submitOption">
                    {{ isSubmited ? `Submitted` : `Submit` }}
                </v-btn>
            </div>
        </template>

    </v-form>
</template>

<script lang="ts" setup>
import { computed, ref, onMounted } from "vue";
import { useChatStore } from '@/stores/chats';
import { type Message, type Operation, type ProtocolDesignInfo, type ProtocolDesignResultInfo, ProtocolDesignFunctionType } from "@/stores/chat-types";
import { uploadSequenceFiles } from '@/api';

// 沟通历史
const chatStore = useChatStore();
// 接口返回的整个对象
const props = defineProps<Message>();
// 接口返回的optionInfo属性
const optionInfo = ref<ProtocolDesignResultInfo>();
// 接口返回的operations属性
const operations = ref<Operation[]>([]);
// 当前步骤
const stage = ref(1);
// 是否已提交
const isSubmited = ref(false);
// 已经上传的文件
const uploadedFiles = ref<string[]>([]);
// 用户上传的文件
const selectedFiles = ref<File[]>([]);
// 文件格式后缀
const fileExtensionRegex = /\.(pdf|txt|doc|docx)$/i;
// 文件格式校验规则
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
// 所有模版的可选项
const protocolSelections = ref<Array<string>>([]);

// 校验用户选择与输入
const optionRules = [
    (value: any) => {
        if (value) return true
        return 'must be select or input parameter'
    },
]


// 文件变化后上传文件
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

// 回调函数, 在组件挂载完成后执行
onMounted(() => {
    optionInfo.value = props?.optionInfo;
    operations.value = props?.optionInfo?.operations || [];
    isSubmited.value = props?.optionInfo?.submitted ? true:false
    stage.value = props?.optionInfo?.stage || 1;
    operations.value.forEach(operation => {
        // file选项
        if (operation.type.includes('file')) {
            uploadedFiles.value = [...operation.value];
        }
        // 加载select项目
        else {
            // 所有模版选项
            protocolSelections.value = [...operation.options];
        }
    });
})

// 是否可以提交
const canSubmit = computed(() => {
    // 文字不在生成中 && 未提交 && (用户选择了协议 || 用户上传了文件)
    const isModified = operations.value.some(operation => operation.value.length > 0);
    return !chatStore.currentChat.isGenerating && !isSubmited.value && isModified;
})

// 提交交互内容
const submitOption = async () => {
    isSubmited.value = true;
    // 当前步骤
    let userReply = "I have submitted";

    let data: ProtocolDesignInfo = {
        stage: props.optionInfo?.stage ?? 2,
        protocol_design_type: props?.optionInfo?.protocol_design_type ?? ProtocolDesignFunctionType.template_protocol_design,
        operations: operations.value,
        // 上一次返回的值
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
