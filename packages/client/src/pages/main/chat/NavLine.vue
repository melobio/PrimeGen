<template>
      <div
      @click="chatStore.hideChatList = !chatStore.hideChatList"
      :class="`nav-line nav-line-${chatStore.hideChatList ? `hide` : `show`
      }`"
    >
      <div class="nav-line-item nav-line-0"></div>
      <div class="nav-line-item nav-line-1"></div>
    </div>
</template>
<script setup lang="ts">
import { useChatStore } from '@/stores/chats';
const chatStore = useChatStore();

</script>
<style lang="scss" scoped>
.nav-line {
    width: 12px;
    height: 30px;
    position: absolute;
    left: 10px;
    top: 50%;
    transform: translateY(-50%);
    cursor: pointer;
    z-index: 100;
    transition: all 0.1s ease;

    &-item {
      width: 5px;
      background-color: #45454e;
      left: 4px;
      position: absolute;
      transition: all 0.1s ease;
      transform-origin: 50% 0;
      height: 18px;
    }
    &-0 {
      top: 0;
      border-top-left-radius: 5px;
      border-top-right-radius: 5px;
    }
    &-1 {
      bottom: 0;
      border-bottom-left-radius: 5px;
      border-bottom-right-radius: 5px;
    }
    &.nav-line-show:hover {
      .nav-line-0 {
        background-color: #666666;
        transform-origin: 50% 0%;
        transform: rotate(17deg) translateY(1px);
        border-bottom-left-radius: 10px;
      }
      .nav-line-1 {
        background-color: #666666;
        transform-origin: 50% 100%;
        transform: rotate(-17deg) translateY(-1px);
        border-top-left-radius: 10px;
      }
    }

    &.nav-line-hide:hover {
      .nav-line-0 {
        transform-origin: 50% 0%;
        transform: rotate(-17deg) translateY(1px);
        border-top-left-radius: 4px;
        border-top-right-radius: 4px;
        border-bottom-right-radius: 10px;
      }
      .nav-line-1 {
        transform-origin: 50% 100%;
        transform: rotate(17deg) translateY(-1px);
        border-bottom-left-radius: 4px;
        border-bottom-right-radius: 4px;
        border-top-right-radius: 10px;
      }
    }
  }
</style>