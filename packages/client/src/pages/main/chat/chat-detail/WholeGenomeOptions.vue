<template>
    <template v-if="operations.length > 0 && operations[0].widget == 'table'">
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
                                <td v-for="head in getTableHeads()" :key="head" class="text-center" style="white-space:nowrap">
                                    <template v-if="head.toLowerCase() == 'ftp path'">
                                        <a :href="item[head.replace(' ', '_').toLowerCase()]" download target="_blank" rel="noopener noreferrer">Detail</a>
                                    </template>
                                    <template v-else >
                                        {{ item[head.replace(/ /g, '_').toLowerCase()] }}
                                    </template>
                                </td>
                            </tr>
                        </tbody>
                    </v-table>
                </div>
            </div>
        </v-form>
    </template>

    <template v-else v-for="operation in operations" :key="operation.key" >

        <template v-if="operation.type.includes('single') && operation.type.includes('select')">
            <div class="option-item">
                <span class="option-item-title">
                    {{ operation.title }}
                </span>
                <span>
                    <v-select size="large" class="option-item-opt-select" v-model="operation.value[0]" clearable
                        :disabled="isSubmited" :items="operation.options">
                    </v-select>
                </span>
            </div>
        </template>

        <template v-if="(operation.type.includes('file'))">
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
                        label="Upload Protocol File" :disabled="isSubmited" >
                    </v-file-input>
                </template>
            </div>
        </template>
    </template>

    <template v-if="props?.optionInfo?.state !== 'stop'">
        <div style="display: flex;justify-content: flex-end;">
            <v-btn variant="elevated" class='ml-3 mt-3' color='#2db489' @click="submitOption"
                :disabled="!canSubmit">
                {{ isSubmited ? `Submited` : `Submit` }}
            </v-btn>
        </div>
    </template>


</template>

<script lang="ts" setup>
import { computed, ref, onMounted, watch } from "vue";
import { useChatStore } from '@/stores/chats';
import { uploadSequenceFiles } from '@/api';
import { type Message, type Operation } from "@/stores/chat-types";

const chatStore = useChatStore();
const props = defineProps<Message>();
const operations = ref<Operation[]>([]);
const uploadedFiles = ref<string[]>([]);
const selectedFiles = ref<File[]>([]);
const isSubmited = ref(false);

const selected = ref<Array<any>>([]);

const selectedChange = (index: number, operation_type: Array<string>) => {
    if (operation_type.includes('single')) {
        selected.value = [index]
    } else {
        selected.value.push(index)
    }
}

const formatColumnName = (name: string) => {
  return name.replace(/(?:^|\s)\w/g, match => match.toUpperCase());
}

const getTableHeads = () => {
    return [
        "accession",
        "representative",
        "organism name",
        "isolate",
        "assembly level",
        "contig n50",
        "scaffold n50",
    ]
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

onMounted(() => {
    operations.value = props?.optionInfo?.operations || [];
    isSubmited.value = props?.optionInfo?.submitted ? true : false;
    operations.value.forEach(operation => {
        if (operation.type.includes('file')) {
            uploadedFiles.value = [...operation.value];
        }
    });
})

const canSubmit = computed(() => {
    const isModified = operations.value.some(operation => operation.value.length > 0);
    return (!chatStore.currentChat.isGenerating) && !isSubmited.value && (selected.value.length > 0 || isModified);
})

const submitOption = async () => {

    isSubmited.value = true;
    let stage = props.optionInfo?.stage ?? 2;
    let userReply = "I have submitted";
    let upload_file_flag = false;

    let search_type = props.optionInfo?.search_type || 'whole_genome_type';
    operations.value.forEach((item, index) => {
        if (!Array.isArray(item.value)) {
            operations.value[index] = { ...item, value: [item.value] }
        }
        if (item.type.includes('file')&&item.value.length>0) {
            upload_file_flag = true;
        }
        else if (selected.value) {
            operations.value[index] = { ...item, value: [item.options[selected.value[0]]] }
        }
    })
    let data = {
        stage,
        search_type,
        operations: operations.value,
        upload_file_flag,
        data: props?.optionInfo?.data ?? {},
    }
    console.log(data);
    await chatStore.sendSelectOptionToNcbiSearch(userReply, data, props)
}

</script>