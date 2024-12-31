<template>
    <v-sheet class="cancer-info-options" style="border-radius: 24px;" width="100%" v-if="operations.length > 0">
        <v-form ref="formRef">
            <template v-for="operation in operations" :key="operation.key">
                <template v-if="(operation.type.includes('single') || operation.type.includes('multi')) && !wantUploadFile">
                    <div class="option-item">
                        <span class="option-item-title">
                            Primary Site:
                        </span>
                        <span class="option-item-opt">
                            <v-select size="large" clearable class="option-item-opt-select" v-model="primarySite"
                                :rules="optionRules" :disabled="isSubmited" :items="primarySiteOptions">
                            </v-select>
                        </span>
                    </div>
                    <div class="option-item">
                        <span class="option-item-title">
                            Site subtype 1:
                        </span>
                        <span class="option-item-opt">
                            <v-select class="option-item-opt-select" v-model="siteSubType1" :rules="optionRules"
                                :disabled="isSubmited" :items="siteSubType1Options"></v-select>
                        </span>
                    </div>
                    <div class="option-item">
                        <span class="option-item-title">
                            Primary histology:
                        </span>
                        <span class="option-item-opt">
                            <v-select class="option-item-opt-select" v-model="primaryHistology" :rules="optionRules"
                                :disabled="isSubmited" :items="primaryHistologyOptions"></v-select>
                        </span>
                    </div>
                    <div class="option-item">
                        <span class="option-item-title">
                            Histology subtype 1:
                        </span>
                        <span class="option-item-opt">
                            <v-select class="option-item-opt-select" v-model="histologySubtype1" :rules="optionRules"
                                :disabled="isSubmited" :items="histologySubtype1Options"></v-select>
                        </span>
                    </div>
                </template>
                <template v-if="operation.type.includes('file') && wantUploadFile">
                    <div class="option-item">
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
                                label="Upload Gene File" :disabled="isSubmited" :rules="fileRules">
                            </v-file-input>
                        </template>
                    </div>
                </template>
            </template>
        </v-form>
        <div class="bottom-btns">
            <v-btn variant='text' v-if="!isSubmited" @click="handleChangeMode">i want
                {{ wantUploadFile ? 'select' : 'upload File' }}</v-btn>
            <v-btn :disabled="!canSubmit" variant="elevated" type="submit" class='ml-3' @click="submitOption">
                {{ isSubmited ? `Submited` : `Submit` }}
            </v-btn>
        </div>
    </v-sheet>
</template>

<script setup lang="ts">
import { watch, onMounted, ref, computed } from 'vue'
import { useChatStore, } from '@/stores/chats';
import { type Message, type CancerOptionInfo, NCBIFunctionType, type Operation } from '@/stores/chat-types';
import { uploadSequenceFiles } from "@/api";
const selectedFiles = ref<File[]>([]);
const uploadedFiles = ref<string[]>([]);
const formRef = ref();
const fileExtensionRegex = /\.(fna|fa)$/i;
const fileRules = ref([
    (files: any) => {
        let validate = true;
        let tips = ''
        files.forEach((file: any) => {
            if (!fileExtensionRegex.test(file.name)) {
                validate = false;
                tips = 'You can only upload fna/fa files'
            }
        });
        if(files.length<1){
            validate = false;
            tips = 'You should upload fna/fa files'
        }
        return validate || tips;
    },
],);
const wantUploadFile = ref<boolean>(false);
const props = defineProps<Message>();
const operations = ref<Operation[]>([]);

const chatStore = useChatStore();
const isSubmited = ref(false);

const cancerOptionInfo = ref()

const primarySiteOptions = ref<Array<string>>([]);
const primarySite = ref();

const siteSubType1Options = ref<Array<string>>([]);
const siteSubType1 = ref();

const primaryHistologyOptions = ref<Array<string>>([]);
const primaryHistology = ref();

const histologySubtype1Options = ref<Array<string>>([]);
const histologySubtype1 = ref();

const optionRules = [
    (value: any) => {
        if (value) return true
        return 'must be selected.'
    },
]

const canSubmit = computed(() => {
    return (!chatStore.currentChat.isGenerating) && (!isSubmited.value);
})

const finishSelected = computed(() => {
    return primarySite.value && siteSubType1.value && primaryHistology.value && histologySubtype1.value ? true : false
})

const finishedUploadFile = computed(() => {
    let fileValidate = false;
    selectedFiles.value.forEach((file: File) => {
        fileValidate = fileExtensionRegex.test(file.name);
    });
    return fileValidate
})

watch(() => primarySite, (newVal, odlVal) => {
    if (newVal.value) {
        siteSubType1Options.value = Object.keys(cancerOptionInfo.value[newVal.value])
    }
    if (!isSubmited.value) {
        siteSubType1.value = '';
        primaryHistology.value = '';
        histologySubtype1.value = '';
    }
}, { deep: true });

watch(() => siteSubType1, (newVal, odlVal) => {
    if (newVal.value) {
        primaryHistologyOptions.value = Object.keys(cancerOptionInfo.value[primarySite.value][newVal.value])
    }
    if (!isSubmited.value) {
        primaryHistology.value = '';
        histologySubtype1.value = '';
    }
}, { deep: true });


watch(() => primaryHistology, (newVal, odlVal) => {
    if (newVal.value) {
        histologySubtype1Options.value = cancerOptionInfo.value[primarySite.value][siteSubType1.value][newVal.value]
    }
    if (!isSubmited.value) {
        histologySubtype1.value = '';
    }
}, { deep: true });

const validate = async () => {
    console.log('formRef.value===>', formRef.value)
    const { valid } = await formRef.value.validate()
    return valid
}


onMounted(() => {
    isSubmited.value = props?.optionInfo?.submitted ? true : false
    if (props?.optionInfo?.operations) {
        operations.value = props?.optionInfo?.operations;
        operations.value.forEach((operation) => {
            if (operation.type.includes('file')) {
                uploadedFiles.value = [...operation.value];
                if (operation.value.length > 0) wantUploadFile.value = true;
            }
            if (operation.type.includes('single') || operation.type.includes('multi')) {
                cancerOptionInfo.value = operation.options || {};
                primarySiteOptions.value = Object.keys(cancerOptionInfo.value);
                primarySite.value = operation.value[0];
                siteSubType1.value = operation.value[1];
                primaryHistology.value = operation.value[2];
                histologySubtype1.value = operation.value[3];
                if (operation.value.length > 0) wantUploadFile.value = false;
            }
        })

    }
})
const handleChangeMode = () => {
    wantUploadFile.value = !wantUploadFile.value;
    selectedFiles.value = [];
    primarySite.value = '';
    siteSubType1.value = '';
    primaryHistology.value = '';
    histologySubtype1.value = '';
    operations.value.forEach((item, index) => {
        operations.value[index].value = [];
    })
}
const handleFileChange = () => {
    if (finishedUploadFile.value) {
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
}
const submitOption = async () => {

    if (await validate()) {
        isSubmited.value = true;
        let stage = props.optionInfo?.stage || 2;
        let userReply = "I have submitted";
        let selected_options = [primarySite.value, siteSubType1.value, primaryHistology.value, histologySubtype1.value]
        operations.value.forEach((item, index) => {
            if (item.key == 'cancer_select') {
                operations.value[index] = { ...item, value: selected_options }
            }
        })
        let cancer_dict = props.optionInfo?.cancer_dict || {};
        let search_type = props.optionInfo?.search_type || NCBIFunctionType.cancer_type;
        let data: CancerOptionInfo = { stage, cancer_dict, search_type, operations: operations.value, data: props?.optionInfo?.data ?? {}, }
        if (finishSelected.value) {
            userReply = ` Primary Site: ${primarySite.value};\n Site subtype 1: ${siteSubType1.value};\n Primary histology: ${primaryHistology.value};\n Histology subtype 1: ${histologySubtype1.value}`;
        }
        await chatStore.sendSelectOptionToNcbiSearch(userReply, data, props)
    }

}
</script>

<style lang="scss">
.cancer-info-options {
    padding: 10px 20px;
    box-sizing: border-box;
    border: 3px solid transparent;
    overflow: hidden;
    margin: 10px 0;

    .option-item {
        width: 100%;
        margin-top: 5px;

        &-title {
            font-size: 14px;
            font-weight: 600;
            margin-bottom: 5px;
            display: block;
        }

        &-opt {

            &-select {}

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

    .bottom-btns {
        margin-top: 10px;
        display: flex;
        justify-content: flex-end;

    }
}
</style>