<template>
    <div class="top-status">
        <div class="status-title">AlphaTool Status:</div>
        <div class="status-content">{{ props.optionInfo?.data?.experiment_info?.status }}</div>

        <!-- commands -->
        <template v-if="canSubmit && props?.optionInfo?.data?.command_options">
            <div class="multi-button-flex">
                <template v-for="command in props?.optionInfo?.data?.command_options" :key="command">
                    <v-btn :disabled="!canSubmit" variant="elevated" class='ml-3 mt-3' color='#2db489' @click="() => submitOption(command)">
                        {{ command }}
                    </v-btn>
                </template>
            </div>
        </template>

    </div>
</template>

<script lang="ts" setup>
import { type Message } from '@/stores/chat-types';
import { useChatStore } from '@/stores/chats';
import { computed, onMounted, ref } from 'vue';

const props = defineProps<Message>();

const isSubmitted = ref(false);

const chatStore = useChatStore();

const canSubmit = computed(() => {
    // not generating
    return !chatStore.currentChat.isGenerating
})

const submitOption = async (command: string) => {
    isSubmitted.value = true;
    let userReply = `${command}`

    await chatStore.sendText(userReply);
}


</script>

<style scoped lang="scss">

.multi-button-flex {
    display: flex;
    justify-content: space-around;
}

.top-status {
    padding: 10px 20px;
    box-sizing: border-box;
    border-radius: 24px;
    border: 3px solid transparent;
    background: #151728;
    margin-top: 20px;

    .status-title {
        color: #afafaf;
        font-family: times, Times New Roman, times-roman, georgia, serif;
        line-height: 30px;
        font-size: 25px;
        text-align: left;
    }

    .status-content{
        color: #ffffff;
        font-family: times, Times New Roman, times-roman, georgia, serif;
        line-height: 45px;
        font-size: 45px;
        font-weight: bold;
        text-align: -webkit-center;
    }

}
</style>
