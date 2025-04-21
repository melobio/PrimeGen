<template>
  <div class="upload-primer">
    <div class="upload-section">
      <div class="file-input-wrapper">
        <input type="file" @change="handleFileSelection" multiple ref="fileInputRef" :disabled="isUploading" style="display: none;" />
      </div>
      <button @click="triggerFileInput" :disabled="isUploading">select primer excel</button>
      <ul>
        <li v-for="(file, index) in selectedFiles" :key="index">{{ file.name }}</li>
      </ul>
      <button @click="uploadFiles" :disabled="selectedFiles.length === 0 || isUploading">upload primer excel</button>
      <div v-if="isUploading" class="loading-animation">
        Loading...
      </div>
    </div>
  </div>
</template>

<script setup lang="ts">
import { deleteMessageById, uploadPrimerExcel } from '@/api';
import { type Message } from '@/stores/chat-types';
import { useChatStore } from '@/stores/chats';
import { reactive, ref } from 'vue';

const props = defineProps<Message>()
const chatStore = useChatStore();
const isUploading = ref(false);

console.log('upload primer: ', props)

const selectedFiles = reactive([] as File[]);
const fileInputRef = ref<HTMLInputElement>();

const triggerFileInput = () => {
  fileInputRef.value?.click();
};

const handleFileSelection = (event: Event) => {

  selectedFiles.splice(0);
  const input = event.target as HTMLInputElement;
  if (input.files) {
    // 将选中的文件追加到selectedFiles数组中
    for (const file of input.files) {
      selectedFiles.push(file);
    }
  }
};

const uploadFiles = async () => {

  isUploading.value = true;

  const formData = new FormData();
  for (const file of selectedFiles) {
    formData.append('files', file);
  }

  try {
    const response = await uploadPrimerExcel(formData);

    if (response.code == 200) {
      const experimentId = response.data as string;

      //save primer and save database
      chatStore.cvstoexper.push({
          conversationUUID: props.conversationUUID,
          experimentId
      });

      let messageEntity = {
        id: props.id,
      }
    
      let messageRet = await deleteMessageById(messageEntity);
      console.log(messageRet);

      //delete cur message
      for (let i = chatStore.currentChatMessages.length - 1; i >= 0; i--) {
        if (chatStore.currentChatMessages[i].id === props.id) {
          chatStore.currentChatMessages.splice(i, 1);
        }
      }

      console.log('upload file start to send text')
      //user say upload ok
      await chatStore.sendText(`Primer info successfully uploaded,experimentId=${experimentId}`);


    }
    
    selectedFiles.splice(0);
    if (fileInputRef.value) {
      fileInputRef.value.value = '';
    }
  } catch (error) {
    console.error('upload file fail', error);
  } finally {
    isUploading.value = false;
  }
};



</script>

<style lang="scss" scoped>
.upload-primer {
  .upload-section {
    background-color: #333;
    padding: 20px;
    border-radius: 10px;
    box-shadow: 0 8px 16px rgba(0, 0, 0, 0.5);
    margin: 10px 0;
    color: #ddd;
    text-align: center;
    transition: background-color 0.3s;

    .file-input-wrapper {
        position: relative;
        display: inline-block;

        &:after {
            color: white;
            background-color: #5cacee;
            padding: 10px 20px;
            border-radius: 20px;
            box-shadow: 0 4px 8px rgba(0, 0, 0, 0.2);
            position: absolute;
            left: 50%;
            top: 50%;
            transform: translate(-50%, -50%);
            transition: background-color 0.3s, box-shadow 0.3s;
            cursor: pointer;
        }

        input[type="file"] {
            position: absolute;
            width: 100%;
            height: 100%;
            left: 0;
            top: 0;
            opacity: 0;
            cursor: pointer;
        }
    }

    button {
        background-color: #48d1cc;
        color: white;
        border: none;
        padding: 10px 20px;
        border-radius: 20px;
        margin-top: 15px;
        cursor: pointer;
        transition: background-color 0.3s, transform 0.3s, box-shadow 0.3s;
        font-weight: bold;
        &:hover {
            background-color: #20b2aa;
            box-shadow: 0 5px 10px rgba(0, 0, 0, 0.3);
            transform: translateY(-2px);
        }
        &:disabled {
            background-color: #aaa;
            cursor: not-allowed;
        }
    }

    ul {
        list-style: none;
        padding: 0;
        margin: 0;
        text-align: left; 
        display: inline-block;
        width: 100%;
        li {
            background-color: #444;
            padding: 5px 10px;
            border-radius: 5px;
            margin-bottom: 5px;
            box-shadow: 0 2px 4px rgba(0, 0, 0, 0.2);
        }
    }
    .loading-animation {
      // 简单的文字加载动画
      &::after {
        content: '.';
        animation: dots 1s steps(5, end) infinite;
      }

      // 定义动画效果
      @keyframes dots {
        0%, 20% {
          color: rgba(0,0,0,0);
          text-shadow:
            .25em 0 0 rgba(0,0,0,0),
            .5em 0 0 rgba(0,0,0,0);
        }
        40% {
          color: black;
          text-shadow:
            .25em 0 0 rgba(0,0,0,0),
            .5em 0 0 rgba(0,0,0,0);
        }
        60% {
          text-shadow:
            .25em 0 0 black,
            .5em 0 0 rgba(0,0,0,0);
        }
        80%, 100% {
          text-shadow:
            .25em 0 0 black,
            .5em 0 0 black;
        }
      }
    }
  }
}
</style>