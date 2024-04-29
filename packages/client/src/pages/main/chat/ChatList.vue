<template>
  <div class='chat-list'>
    <div class='search-container'>
      <div class='search'>
        <v-icon icon='mdi-magnify' style='color: #545665'/>
        <input class='search-input' placeholder='Search for Chats'/>
      </div>
      <v-btn icon="mdi-plus" size='40px' class='ml-3' rounded color='#1E2032' @click="addNewConWithExp"></v-btn>
    </div>
    <div class='list-container'>
      <chat-list-item
        v-for='(chat) in chatStore.chats'
        v-bind='chat'
        @click='handleSelectedChat(chat)'
        :key='chat.uuid'
      />
    </div>
  </div>
</template>

<script setup lang='ts'>
import ChatListItem from '@/pages/main/chat/ChatListItem.vue';
import type { Chat } from '@/stores/chat-types';
import { useChatStore } from '@/stores/chats';
import { useExperimentsStore } from '@/stores/experiments';
const chatStore = useChatStore();
const experimentsStore = useExperimentsStore();

const addNewConWithExp =async () =>{
  const res =await experimentsStore.handleAddConWithExp()
  if(res.code == 200){
    await experimentsStore.getAllExperiments();
    await chatStore.getAllChats();
    chatStore.currentChat = chatStore.chats[chatStore.chats.length-1];
  }
}
const handleSelectedChat = (chat:Chat)  => {
  chatStore.currentChat = chat
}
</script>

<style scoped lang='scss'>
.chat-list {
  width: 250px;
  min-width: 250px;
  height: 100%;
  padding: 15px 0;
  display: flex;
  flex-direction: column;
  .search-container {
    //padding: 0 15px;
    padding-left: 15px;
    margin-top: 25px;
    display: flex;
    flex-direction: row;
    .search {
      flex: 1;
      height: 40px;
      background: #1E2032;
      border-radius: 6px;
      display: flex;
      align-items: center;
      padding: 0 15px;
      .search-input {
        height: 100%;
        flex: 1;
        margin-left: 10px;
        border: none;
        width: 0px;
        outline: none;
        color: white;
        &::placeholder {
          color: #545665;
        }
      }
    }
  }

  .list-container {
    margin-top:10px;
    padding: 20px 0 0 15px;
    overflow-y: auto;
    flex: 1;
  }
}
</style>