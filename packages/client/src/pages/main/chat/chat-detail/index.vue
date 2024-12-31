<template>
  <div class='chat-detail-container'>
    <div class='chat-title'>
      <div class="chat-title-text">Experiment of {{chatStore.currentChat?.name}} </div>
      <img v-if="!isFullScreen" style='width: 25px; height: 25px;' @click="handleFullscreen" src='@/assets/fullscreen.svg'
        alt='' />
      <!-- <el-button plain style='width: 25px; height: 25px;margin-left:20px' @click="displayVideos" >
        video
      </el-button> -->
    </div>
    <div class='chat-list' ref="chatLists" @scroll="scroll">
      <template v-for="(message) in chatStore.currentChatMessages" :key="message.id">
        <Remote v-bind="message" v-if="message.role === 'assistant'" />
        <Local v-bind="message" v-else-if="message.role === 'user'" />
      </template>
    </div>
    <div class='chat-input'>
      <textarea class='input' v-model.trim="inputText" placeholder="Message PrimeGen..." @keydown="onInputKeyDown"
        autofocus ref="textareaInput" />
      <v-btn @click="stopGenerating = true;" v-if="!canSend" color="transparent" size="samll"
        style="position: absolute; top: 50%; right: 180px; transform: translateY(-50%);">
        <v-progress-circular indeterminate color="#31DA9F">
          <v-icon color="#31DA9F" icon="mdi mdi-stop" size="small" />
        </v-progress-circular>
      </v-btn>
      <AudioRecordButton :disabled="!canSend" />
      <v-btn @click="send" elevation='0' :disabled="!canSend" size='60' color='#292B3C' style='border-radius: 10px;'>
        <img src='@/assets/send.png' style='width: 25px; height: 25px;' alt='' />
      </v-btn>
    </div>
    <TipsDialog />
    <RestartDialog />
    <StopGeneratingDialog v-model:stopGenerating="stopGenerating" />
    <NavLine />
    <!-- <div class="videoContainer" v-if="dialogVisible" >
      <VideoDisplay/>
    </div> -->
  </div>
</template>

<script setup lang='ts'>
import { useChatStore } from '@/stores/chats';
import NavLine from '@/pages/main/chat/NavLine.vue'
import Remote from '@/pages/main/chat/chat-detail/Remote.vue'
import Local from '@/pages/main/chat/chat-detail/Local.vue'
import { nextTick, onMounted, ref, watch, computed } from "vue";
import AudioRecordButton from './AudioRecordButton.vue'
import TipsDialog from '@/pages/main/chat/chat-detail/TipsDialog.vue'
import RestartDialog from '@/pages/main/chat/chat-detail/RestartDialog.vue'
import StopGeneratingDialog from '@/pages/main/chat/chat-detail/StopGeneratingDialog.vue'
// import VideoDisplay from '@/pages/main/chat/chat-detail/VideoDisplay.vue'
// import { ElMessageBox } from 'element-plus'
const chatStore = useChatStore();
const inputText = ref('');
const chatLists = ref();
const textareaInput = ref<HTMLButtonElement | null>(null);
const stopGenerating = ref<boolean>(false);

const isFullScreen = ref(false)
onMounted(() => {
  chatStore.scrollBottom = true;
});

const canSend = computed(() => {
  let canSend = true;
  let hasGeneratingMsg = chatStore.currentChatMessages.filter(item => item.generating)
  if (chatStore.isRecordingAudio || chatStore.currentChat.isGenerating || hasGeneratingMsg.length > 0) {
    canSend = false;
  }
  return canSend
})

watch(() => chatStore.scrollBottom, (scrollBottom) => {
  if (scrollBottom) {
    nextTick(() => {
      chatLists.value.scrollTo(0, chatLists.value.scrollHeight)
    })
    chatStore.scrollBottom = false
  }
})

watch(() => chatStore.currentChat.id, (newValue, oldValue) => {
  if (textareaInput.value) {
    textareaInput.value.focus();
  }
});
const dialogVisible = ref(false);
const displayVideos = () => {
dialogVisible.value = !(dialogVisible.value);
console.log(dialogVisible.value);
}

const handleFullscreen = () => {
  if (!document.fullscreenElement) {
    document.documentElement.requestFullscreen();
    isFullScreen.value = true;
    chatStore.hideChatList = true;
  } else {
    if (document.exitFullscreen) {
      document.exitFullscreen();
      isFullScreen.value = false;
      chatStore.hideChatList = false;
    }
  }
}
async function send(ev: KeyboardEvent) {
  if (!inputText.value) {
    console.warn('input text is empty');
    return;
  }
  if (!canSend.value) {
    console.warn('can`t Send');
    return
  }
  if (!ev.shiftKey) {
    await chatStore.sendText(inputText.value);
    ev.preventDefault();
    inputText.value = '';
  }
}
function onInputKeyDown(ev: KeyboardEvent) {
  if (ev.key === "Enter") {
    ev.preventDefault();
    if (ev.ctrlKey) {
      inputText.value += "\n";
    } else {
      send(ev)
    }
  }
}

function scroll(ev: any) {
  chatStore.shouldScrollToBottomWhileGenerating =
    (ev.target.scrollHeight - ev.target.clientHeight) - ev.target.scrollTop <= 1;
}
</script>

<style scoped lang='scss'>
.chat-detail-container {
  background-color: #1E2032;
  width: 100%;
  flex: 1;
  display: flex;
  flex-direction: column;
  overflow: hidden;
  position: relative;

  .videoContainer {
    position: absolute;
    top: 0%;
    left: 50%;
    transform: translateX(-50%)
  }

  .chat-title {
    height: 60px;
    max-width: 750px;
    line-height: 60px;
    font-size: 16px;
    padding: 0px 20px;
    text-overflow: ellipsis;
    overflow: hidden;
    display: flex;
    align-items: center;
    white-space: nowrap;
    color: #ECECEC;
    border-bottom: 1px solid #292B3C;

    &-text {
      max-width: 90%;
      overflow: hidden;
      text-overflow: ellipsis;
    }

    img {
      margin-left: 10px;
    }
  }

  .chat-list {
    flex: 1;
    height: 100%;
    overflow-y: auto;
  }

  .chat-input {
    margin-top: 20px;
    border-top: 1px solid #292B3C;
    position: relative;
    padding: 20px;
    box-sizing: content-box;
    display: flex;
    align-items: stretch;

    .input {
      background-color: #151728;
      border-radius: 10px;
      flex: 1;
      outline: none;
      padding: 10px;
      color: #ECECEC;
      resize: none;
      line-height: 20px;
      margin-right: 15px;

      &[disabled] {
        cursor: not-allowed;
      }
    }
  }
}
</style>