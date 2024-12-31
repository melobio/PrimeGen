<template>
  <div class='chat-list' :style="{maxWidth:chatStore.hideChatList?'0px':'250px'}">
    <div class='search-container'>
      <div class='search'>
        <v-icon icon='mdi-magnify' style='color: #545665'/>
        <input class='search-input' v-model="searchInput" placeholder='Search for Chats'/>
        <span class="close-icon" v-if="searchInput" @click="searchInput = ''">
          <v-icon icon='mdi mdi-close-circle'/>
        </span>
      </div>
      <v-btn icon="mdi-plus" size='40px' class='ml-3' rounded color='#1E2032' @click="addNewConWithExp"></v-btn>
    </div>
    <div class='list-container' ref="chatListContainer">
      <chat-list-item
        v-for='(chat) in chatList'
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
import {nextTick,watch,ref,computed} from "vue";
const chatListContainer = ref();
const searchInput = ref('');
const chatStore = useChatStore();
const experimentsStore = useExperimentsStore();
const chatList = computed(()=>{
  return chatStore.chats.filter(item => item.name.includes(searchInput.value))
})
const addNewConWithExp =async () =>{
  const res =await experimentsStore.handleAddConWithExp()
  if(res.code == 200){
    await experimentsStore.getAllExperiments();
    await chatStore.getAllChats();
    chatStore.currentChat = chatStore.chats[0];
  }
}
watch(() => chatStore.currentChat, (currentChat) => {
  if (currentChat?.id == chatStore.chats[0]?.id) {
    nextTick(() => {
      chatListContainer.value.scrollTop = 0;
    })
  }
})

const handleSelectedChat = (chat:Chat)  => {
  chatStore.currentChat = chat
}
</script>

<style scoped lang='scss'>
.chat-list {
  max-width: 250px;
  height: 100%;
  padding: 15px 0;
  transition: all 0.2s ease-out;
  display: flex;
  flex-direction: column;
  overflow: hidden;
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
      padding: 0 10px;
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
      .close-icon{
        width: 22px;
        margin-left: 10px;
        cursor: pointer;
        
        .mdi-close-circle {
          color: #545665;
          }
        &:hover{
          .mdi-close-circle {
            color:#98989e;
          }
        }
      }
    }
  }

  .list-container {
    margin-top:10px;
    padding: 20px 0 0 15px;
    overflow-y: auto;
    overflow-x: hidden;
    flex: 1;
  }
}
</style>