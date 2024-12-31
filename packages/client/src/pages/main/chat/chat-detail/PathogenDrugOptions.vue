<template>
    <div class="pathogen-drug-options" v-if="operations.length">
        <div>
            <v-table style="background: transparent;">
                <thead>
                    <tr>
                        <th class="text-center">
                            Gene
                        </th>
                        <th class="text-center">
                            Position
                        </th>
                        <th class="text-center">
                            Drug
                        </th>
                        <th class="text-center">
                            Select
                        </th>
                    </tr>
                </thead>
                <tbody>
                    <template v-if="operations.length">
                        <template v-for="(option) in operations">
                            <tr class="text-center" v-for="(item,index) in option.options" :key="index">
                                <td class="text-center">{{ item.gene }}</td>
                                <td class="text-center">{{ item.position }}</td>
                                <td class="text-center">{{ item.Drug }}</td>
                                <td class="text-center">
                                    <v-checkbox hide-details 
                                        :disabled="isSubmited" 
                                        v-model="selected"
                                        :value="item.Drug"></v-checkbox>
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
<script lang="ts" setup>
import { computed, ref, onMounted } from "vue";
import { useChatStore } from '@/stores/chats';
import { NCBIFunctionType, type Message, type PathogenDrugInfo ,type Operation} from "@/stores/chat-types";
import { AgentType } from "@xpcr/common/src/agent";
const operations = ref<Operation[]>([]);
const chatStore = useChatStore();
const props = defineProps<Message>();
const selected = ref<Array<any>>([]);
const isSubmited = ref(false);

onMounted(() => {
    operations.value = props?.optionInfo?.operations || [];
    isSubmited.value = props?.optionInfo?.submitted ? true:false

    operations.value.forEach(operation => {
        if (!operation.type.includes('file')) {
            operation.options.forEach((item:Record<string,any>) => {
                operation.value.forEach(operationItem => {
                    if (operationItem == item['Drug']) {
                        selected.value.push(operationItem)
                    }
                })
            })
        }
    });
})

const canSubmit = computed(() => {
    // const lastSearchMsgId = chatStore.currentChatMessages.filter(item => item.agentType == AgentType.SEQUENCE_SEARCH).pop()?.id || '';
    // props.id == lastSearchMsgId && 
    return (!chatStore.currentChat.isGenerating) && !(isSubmited.value);
})


const submitOption = async () => {
    isSubmited.value = true;
    let stage = props.optionInfo?.stage ?? 2;
    let pathogen_drug_resistance_dict = props.optionInfo?.pathogen_drug_resistance_dict ?? {};
    let userReply = "I want to study drug" + "'" + selected.value.join("'and drug '") + "'";
    let search_type = props.optionInfo?.search_type || NCBIFunctionType.species_identification_type;
    operations.value.forEach((item,index)=>{
          operations.value[index] = {...item,value:[...selected.value]}
    })
    let data: PathogenDrugInfo = {
        stage,
        search_type,
        pathogen_drug_resistance_dict,
        operations: operations.value,
        // last response data
        data: props?.optionInfo?.data ?? {},
    }
    await chatStore.sendSelectOptionToNcbiSearch(userReply, data, props)
}
</script>
<style scoped lang="scss">
.pathogen-drug-options {
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
