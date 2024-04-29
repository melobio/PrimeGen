<template>
    <div class="initiative-start-primer-design" v-if="canSubmit">
        <div style="display: flex;justify-content: flex-end;">
            <v-btn :disabled="chatStore.currentChat.isGenerating" variant="elevated" class='ml-3 mt-3' color='#2db489' @click="submit">
                start primer design
            </v-btn>
        </div>
    </div>
</template>
<script lang="ts" setup>
import { computed } from "vue";
import { useChatStore } from '@/stores/chats';
import { type Message } from "@/stores/chat-types";

const chatStore = useChatStore();
const props = defineProps<Message>();

const canSubmit = computed(() => {
    const lastSearchMsgId = chatStore.currentChatMessages.slice().pop()?.id || '';
    return  props.id === lastSearchMsgId;
})

const submit = async () => {
   chatStore.sendText('please start Primer Design');
}

</script>
<style scoped lang="scss">
.initiative-start-primer-design {
    padding: 10px 20px;
    box-sizing: border-box;
}
</style>
