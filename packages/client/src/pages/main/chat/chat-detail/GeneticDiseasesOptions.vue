
<template>
        <div class="genetic-diseases-options">
          <div>
            <v-table style="background: transparent;">
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
                <tbody >
                    <template  v-if="showOptions">
                        <tr v-for="item in optionList" :key="item.index">
                            <td class="text-center">{{ item.gene }}</td>
                            <td class="text-center">{{ item.disease }}</td>
                            <td class="text-center">
                                <v-checkbox
                                    hide-details
                                    :disabled="isSubmited"
                                    v-model="selected"
                                    :value="item.disease"
                                ></v-checkbox>
                            </td>
                        </tr>
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
           </div>
          <div style="display: flex;justify-content: flex-end;">
            <v-btn :disabled="!canSubmit" variant="elevated" class='ml-3 mt-3' color='#2db489' @click="submitOption">
              {{ isSubmited ? `Submited` : `Submit` }}
            </v-btn>
          </div>
        </div>
  </template>
  
<script setup lang="ts">
import {computed,ref,watch,onMounted} from "vue";
import { useChatStore } from '@/stores/chats';
import type { Message,GeneticDisorderInfo } from "@/stores/chat-types";
import { NCBIFunctionType } from  "@/stores/chat-types";

const chatStore = useChatStore();
const selected = ref<Array<any>>([]);
const props = defineProps<Message>()
const optionList = ref<Array<Record<string,string>>>([])
const isSubmited = ref(false);


const showOptions= computed(()=>{
    return (props?.optionInfo?.options||[]).length>0;
})

const  canSubmit= computed(()=>{
    // const lastSearchMsgId = chatStore.currentChatMessages.filter(item=>item.agentType == AgentType.SEQUENCE_SEARCH).pop()?.id || '';
    // props.id == lastSearchMsgId &&
    return (!chatStore.currentChat.isGenerating) &&  selected.value?.length>0 && !(isSubmited.value);
})

watch(() => props.optionInfo, (newVal, oldVal) => {
  let optionData = newVal?.options || [];
  let options = optionData.map(obj => ({ gene: obj.gene, disease: obj.disease }));
  optionList.value = options;
}, { deep: true, immediate: true });

onMounted(()=>{
  let selectedOptionsInfo = props.optionInfo?.selectedOptionsInfo||[]
  isSubmited.value = selectedOptionsInfo.length>0;
  if(isSubmited.value){
    selected.value = selectedOptionsInfo.map(item=>item.disease)
  }
})

const submitOption = async () => {
  let stage = props.optionInfo?.stage ?? 2;
  let state = props.optionInfo?.state ?? 'continue';
  let search_type =  props.optionInfo?.search_type ?? NCBIFunctionType.genetic_disorder_type;
  let selectedOptionsInfo = selected.value.map((a: string) => {
      const opt = optionList.value.find(b => a === b.disease);
      return opt ? { ...opt } : {};
  });
  let selectedGene = selectedOptionsInfo.map((item:any)=>item.disease).map((a: string) => {
    return optionList.value.find(b => a === b.disease)?.gene || '';
  })
  let selected_gene_options = selectedGene.filter((item: string, index: any) => selectedGene.indexOf(item) === index && item !== '');
  // I want to perform an analysis on the genes 'CHRNA1' and 'CHRND'.
  let options = props?.optionInfo?.options || [];
  let userReply = "I want to perform an analysis on the genes " + "'" + selected_gene_options.join("' and '") + "'";
  let data :GeneticDisorderInfo= { options,stage,state,selected_options:selected_gene_options, search_type}
  isSubmited.value = true;
  await chatStore.sendSelectOptionToNcbiSearch(userReply,data,props)
}
</script>

<style scoped lang="scss">
.genetic-diseases-options {
    padding: 10px 20px;
    box-sizing: border-box;
    border-radius: 24px;
    border: 3px solid transparent;
    background: #151728;
    margin: 10px 0;
    .v-input{
      display: flex;
      justify-content: center;
    }
}
</style>