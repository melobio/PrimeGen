
<template>
  <div class="genetic-diseases-options" v-if="operations?.length > 0">
      <template v-if="optionInfo?.data?.max_length && optionInfo?.data.max_length > 1">
        <v-form ref="formRef">
            <div class="option-item" v-for="(item, index) in operations" :key="item.title">
                <span class="option-item-title">
                    <span class="required-icon" v-if="item?.type.includes('required')">*</span>
                    {{
                        item.key.split('_').join(' ').toUpperCase() }}
                </span>
                <div class="option-item-opt">
                    <v-text-field type="number" v-if="item.type.includes('input')" clearable
                        class="option-item-opt-input" v-model="operations[index].value[0] " :placeholder="`please input`"
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
    <div v-else v-for="(operation, index) in operations" :key="index" class="option-item">
      <v-table style="background: transparent;" v-if="operation.key == 'gene_select'">
        <thead>
          <tr>
            <th class="text-center">
              Gene
            </th>
            <th class="text-center">
              disease
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
                <td class="text-center">{{ item.gene }}</td>
                <td class="text-center">{{ item.disease }}</td>
                <td class="text-center">
                  <v-checkbox hide-details :disabled="isSubmited" v-model="selected" :value="item.gene"></v-checkbox>
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
          </tr>
        </tbody>
      </v-table>
      <template v-if="operation.key == 'gene_path'">
        <div class="option-item" v-for="(item, index) in operations" :key="item.title">
          <span class="option-item-title">
            <span class="required-icon" v-if="item?.type.includes('required')">*</span> {{ item.key.split('_').join('').toUpperCase() }}
          </span>
          <div class="option-item-des">
            {{ item.title }} :
          </div>
          <div class="option-item-opt">
            <template v-if="item.type.includes('file')">
              <template v-if="item.value.length > 0 && isSubmited">
                <div class="file-uploaded" v-for="(fileItem, fileIndex) in operations[index].value" :key="fileIndex">
                  {{ fileItem.split('/').pop() }}
                </div>
              </template>
              <template v-else>
                <v-file-input class="file-input" v-model="selectedFiles[item.key]" @change="handleFileChange(item.key)"
                  :label="`please upload`" prepend-icon="" :disabled="isSubmited" :rules="rules">
                </v-file-input>
              </template>
            </template>
          </div>
        </div>

      </template>
      <!--  -->
      <div style="display: flex;justify-content: flex-end;">
        <v-btn :disabled="!canSubmit" variant="elevated" class='ml-3 mt-3' color='#2db489' @click="submitOption">
          {{ isSubmited ? `Submited` : `Submit` }}
        </v-btn>
      </div>
    </div>
  </div>
</template>
  
<script setup lang="ts">
import { computed, ref, watch, onMounted } from "vue";
import { useChatStore } from '@/stores/chats';
import type { Message, GeneticDisorderInfo, Operation, GeneticLengthInfo } from "@/stores/chat-types";
import { NCBIFunctionType } from "@/stores/chat-types";
import { uploadSequenceFiles } from "@/api";
const geneFileregex = /^(.+?)_gene\.(fna|fa)$/;

const validateRules = [
    {
        keys: ['gene_path'],
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

// 1. start 开始序列
// 2. end 结束序列
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
                let max_length = props.optionInfo?.max_length || Infinity;
                if (value == '' || value > max_length) {
                    return 'End position larger than ' + max_length
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

const selectedFiles = ref<Record<string, any>>({});

const chatStore = useChatStore();
const selected = ref<Array<any>>([]);
const props = defineProps<Message>()
const operations = ref<Operation[]>([]);
const isSubmited = ref(false);

const canSubmit = computed(() => {
  return (!chatStore.currentChat.isGenerating) && selected.value?.length > 0 && !(isSubmited.value);
})

const lengthCanSubmit = computed(() => {
    return (!chatStore.currentChat.isGenerating) && !isSubmited.value;
})

onMounted(() => {
  operations.value = props?.optionInfo?.operations || [];
  isSubmited.value = props?.optionInfo?.submitted ? true : false

  if (isSubmited.value) {
    operations.value.forEach(item => {
      if (item.key == 'gene_select') {
        selected.value = item.value;
      }
    })
  }
})
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
  isSubmited.value = true;
  let stage = props.optionInfo?.stage ?? 2;
  let state = props.optionInfo?.state ?? 'continue';
  let search_type = props.optionInfo?.search_type ?? NCBIFunctionType.genetic_disorder_type;
  let selected_gene_options: any[] = [];

  if (props.optionInfo?.data?.max_length && props.optionInfo?.data?.max_length > 1) {
      let userReply = "I have submitted";
      let start_pos = Number(getOperationValue("start_pos"));
      let end_pos = Number(getOperationValue("end_pos"));
      let target_gene_path = props.optionInfo?.target_gene_path || [];
      let data: GeneticLengthInfo = {
          stage,
          search_type,
          start_pos,
          end_pos,
          target_gene_path,
          // 上一次返回的值
          data: props?.optionInfo?.data ?? {},
      }
      await chatStore.sendSelectOptionToNcbiSearch(userReply, data, props)
  } else {
    operations.value.forEach((operation, index) => {
      if (operation.key == 'gene_select') {
        operations.value[index].value = selected.value
        selected_gene_options = operation.options.filter((item: { disease: any; gene: any}) => selected.value.includes(item.gene)).map((item: { disease: any; }) => item.disease)
      }
    })
    // I want to perform an analysis on the genes 'CHRNA1' and 'CHRND'.
    let userReply = "I want to perform an analysis on the genes " + "'" + selected_gene_options.join("' and '") + "'";
    let data: GeneticDisorderInfo = { operations: operations.value, stage, state, search_type }
    await chatStore.sendSelectOptionToNcbiSearch(userReply, data, props)
  }

}
</script>

<style  lang="scss">
.genetic-diseases-options {
  padding: 10px 20px;
  box-sizing: border-box;
  border-radius: 24px;
  border: 3px solid transparent;
  background: #151728;
  margin: 10px 0;


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
    width: 100%;
    .v-messages {
      .v-messages__message {
        padding: 5px 0;
      }
    }
  }
  .v-input__control{
    justify-content: center;
  }
  .v-selection-control {
    flex: none !important;
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
}
</style>