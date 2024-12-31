<template>
    <v-dialog v-model="tipsDialog" width="auto">
        <div class='restart-dialog'>
            <div class='title'>
                Restart  Experiment
            </div>
            <v-divider :thickness="2" color="warning"></v-divider>
            <div class='content'>
                <div class="content-text">
                    Are you sure you want to restart the Experiment?
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

const chatStore = useChatStore();

const tipsDialog = computed(() => {
    return chatStore.restartConversationUUID !== '';
});

async function handleOperate(restart: boolean) {
    chatStore.handleRestartConversation(restart)
}
</script>
  
<style scoped lang='scss'>
.restart-dialog {
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