<template>

    <template v-for="operation in operations" :key="operation.key" >
        <template v-if="(operation.type.includes('file'))">
            <div class="option-item">
                <span class="option-item-title">
                    {{ operation.title }}
                </span>
                <template v-if="uploadedFiles.length > 0 && isSubmitted">
                    <div class="file-uploaded" v-for="(fileItem, fileIndex) in uploadedFiles" :key="fileIndex">
                        {{ fileItem.split('/').pop() }}
                    </div>
                </template>
                <template v-else>
                    <v-file-input class="file-input" v-model="selectedFiles[operation.key]" @change="handleFileChange(operation.key)"
                        :label="`please upload`" prepend-icon="" :disabled="isSubmitted">
                    </v-file-input>
                </template>
            </div>
        </template>
    </template>

    <template v-if="props?.optionInfo?.state !== 'stop'">
        <div style="display: flex;justify-content: flex-end;">
            <v-btn variant="elevated" class='ml-3 mt-3' color='#2db489' @click="submitOption"
                :disabled="!canSubmit">
                {{ isSubmitted ? `Submited` : `Submit` }}
            </v-btn>
        </div>
    </template>

</template>

<script setup lang="ts">
import { computed, ref, onMounted } from "vue";
import { useChatStore } from '@/stores/chats';
import { uploadSequenceFiles } from '@/api';
import { type Message, type Operation, PrimerDesignFunctionType } from "@/stores/chat-types";

const chatStore = useChatStore();
const props = defineProps<Message>();
const operations = ref<Operation[]>([]);
const uploadedFiles = ref<string[]>([]);
const selectedFiles = ref<Record<string, any>>({});
const isSubmitted = ref(false);


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
    console.log('selectedFiles', selectedFiles.value);
    console.log('operations', operations.value);
}

onMounted(() => {
    operations.value = props?.optionInfo?.operations || [];
    isSubmitted.value = props?.optionInfo?.submitted ? true : false;
    operations.value.forEach(operation => {
        if (operation.type.includes('file')) {
            uploadedFiles.value = [...operation.value];
        }
    });
})

const canSubmit = computed(() => {
    const isModified = operations.value.every(operation => operation.value.length > 0);
    console.log('isModified', isModified);
    return (!chatStore.currentChat.isGenerating) && !isSubmitted.value && isModified;
})

const submitOption = async () => {

    isSubmitted.value = true;
    let stage = props.optionInfo?.stage ?? 2;
    let userReply = "I have submitted";
    let upload_file_flag = false;

    let primer_type: PrimerDesignFunctionType = props.optionInfo?.primer_type || PrimerDesignFunctionType.redesign_primer_type;
    operations.value.forEach((item, index) => {
        if (!Array.isArray(item.value)) {
            operations.value[index] = { ...item, value: [item.value] }
        }
        if (item.type.includes('file')&&item.value.length>0) {
            upload_file_flag = true;
        }
    })
    let data = {
        stage,
        primer_type,
        operations: operations.value,
        upload_file_flag,
        data: props?.optionInfo?.data ?? {},
    }
    await chatStore.sendSelectOptionToPrimerDesign(userReply, data, props)
}

</script>

<style scoped lang="scss">
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
</style>