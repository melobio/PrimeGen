<template>
  <div class='main'>
    <div class='menu-container' :style="{width:chatStore.hideChatList ? '0px':'200px'}" >
      <Menu v-model='menu'/>
    </div>
    <div class='content-container'>
      <Chat v-if="menu === 'Chats'"/>
      <Experiment v-else-if="menu === 'Experiment'"/>
      <Tools v-else-if="menu === 'Tools'" />
    </div>
  </div>
</template>

<script setup lang='ts'>
import Menu from '@/pages/main/Menu.vue'
import Chat from '@/pages/main/chat/index.vue'
import Experiment from "@/pages/main/experiment/index.vue"
import Tools from "@/pages/main/tools/index.vue"
import {onMounted, ref} from 'vue'
import { useChatStore} from "@/stores/chats";
const chatStore = useChatStore();

onMounted(() => {
  chatStore.init();
});

const menu = ref('Chats')
</script>

<style scoped lang='scss'>
.main {
  width: 100%;
  height: 100%;
  display: flex;
  flex-direction: row;
  align-items: stretch;
  background: #151728;
  .menu-container {
    transition: all 0.2s ease-out;
    width: 200px;
  }
  .content-container {
    flex: 1;
  }
}
</style>