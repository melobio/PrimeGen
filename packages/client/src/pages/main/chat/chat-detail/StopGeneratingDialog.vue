<template>
    <v-dialog v-model="stopGenerating" width="auto">
        <div class='stop-generating'>
            <div class='title'>
                Stop  Generating
            </div>
            <v-divider :thickness="2" color="warning"></v-divider>
            <div class='content'>
                <div class="content-text">
                    Are you sure you want to stop  Generating message?
                </div>
            </div>
            <div class='text-right mt-4'>
                <v-btn variant="text" color='#2db489' @click='handleOperate(false)'>
                    <span style='text-transform: none'>
                        Cancel
                    </span>
                </v-btn>
                <v-btn variant='elevated' class='ml-3' color='red' @click='handleOperate(true)'>
                    <span style='text-transform: none'>
                        Confirm
                    </span>
                </v-btn>
                
            </div>
        </div>
    </v-dialog>
</template>
  
<script setup lang='ts'>
import { computed } from "vue";
import { useChatStore } from '@/stores/chats';
const props = defineProps<{
    stopGenerating:boolean
}>()
const emit = defineEmits(['update:stopGenerating'])
const stopGenerating = computed(()=>{
    return props.stopGenerating
})  
const chatStore = useChatStore();


async function handleOperate(stop: boolean) {
    if(stop){
        chatStore.stopGeneratingMessage()
    }
    emit('update:stopGenerating', false)
}
</script>
  
<style scoped lang='scss'>
.stop-generating {
    background-color: #1E2032;
    border-radius: 10px;
    padding: 10px 20px 10px;

    .title {
        color: #ECECEC;
        font-size: 18px;
        line-height: 44px;
        height: 44px;
    }

    .content {
        color: #ECECEC;
        font-size: 15px;
        display: flex;
        justify-content: center;
        align-items: center;
        flex-direction: column;

        .content-text {
            width: 500px;
            line-height: 24px;
            margin-bottom: 10px;
        }
    }
}
</style>