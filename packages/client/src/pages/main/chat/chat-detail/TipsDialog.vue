<template>
  <v-dialog
      v-model="tipsDialog"
      width="auto"
    >
    <div class='tips-dialog'>
    <div class='title'>
       Warning
    </div>
    <v-divider :thickness="2" color="warning"></v-divider>
    <div class='content'>
        <div class="content-text">
          It seems that a bent tip has appeared. The task execution is currently paused. Please check the device, and then you can choose to resume or stop the task.      
        </div>
        <JetFaultCheckResult
            v-bind="chatStore.badTipInfo"
            :img-width="450"
          />
    </div>
    <div class='text-right mt-4'>
      <v-btn variant='text'  color='red' @click='handleOperate("stop")'>
        <span style='text-transform: none'>
          Stop
        </span>
      </v-btn>
      <v-btn variant="elevated" class='ml-3' color='#2db489' @click='handleOperate("play")'>
        <span style='text-transform: none'>
          Resume
        </span>
      </v-btn>
    </div>
  </div>
    </v-dialog>
</template>

<script setup lang='ts'>
import { computed } from "vue";
import { useChatStore } from '@/stores/chats';
import JetFaultCheckResult from "@/pages/main/chat/chat-detail/JetFaultCheckResult.vue";

const chatStore = useChatStore();

const tipsDialog = computed(() => {
  return chatStore.badTipInfo.bad_tip>0;
});

async function handleOperate(optionType:'stop'|'play'){
  let res:{success:boolean,data:string}|undefined= await chatStore.sendOption(optionType,chatStore.badTipInfo.runId);
  if(res?.success){
    chatStore.setBadTipInfo({
    runId:'',
    bad_tip:0,
    pcr:0,
    tip:0, 
    conversationUUID:'' ,
    data:'',
    plot:[]
  })
  }
}
</script>

<style scoped lang='scss'>
.tips-dialog {
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
    .content-text{
      width: 500px;
      line-height: 24px;
      margin-bottom: 10px;
    }
  }
}
</style>