<template>
  <div class='chat-container'>
    <ChatList />
    <ChatDetail
      style='border-top-left-radius: 30px;'
      class='flex-1-1 w-100 mt-10 ml-4'/>
  </div>
</template>

<script setup lang='ts'>
import { onMounted } from 'vue'
import { useChatStore } from '@/stores/chats'
import ChatList from '@/pages/main/chat/ChatList.vue'
import ChatDetail from '@/pages/main/chat/chat-detail/index.vue';

onMounted(() => {
  useChatStore().getAllChats();
});
</script>

<style scoped lang='scss'>
.chat-container {
  display: flex;
  flex-direction: row;
  align-items: stretch;
  width: 100%;
  height: 100%;
}
</style>