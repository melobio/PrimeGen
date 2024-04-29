<template>
  <div class="upload-sequence-file">
    {{ props.content }}
    <div class="upload-section">
      <div style="padding:10px">
        <v-file-input
        v-model="selectedFiles"
        label="Select gene files"
        prepend-icon="mdi-file-upload"
        :rules="rules"
        multiple
        :disabled="isUploading || !canSelect"></v-file-input>
        <div style="text-align: right;">
          <v-btn :disabled="!canUpload" variant="elevated" class='ml-3 loading-animation upload-btn' @click="handleUploadFiles"
          v-if="selectedFiles.length > 0">
          {{isUploading? 'Uploading':'Upload'}}
        </v-btn>  
        </div>
      </div>
    </div>
  </div>
</template>

<script setup lang="ts">
import { uploadSequenceFiles } from '@/api';
import { type Message } from '@/stores/chat-types';
import { useChatStore } from '@/stores/chats';
import { ref ,computed} from 'vue';
import { AgentType } from "@xpcr/common/src/agent";

const props = defineProps<Message>()
const chatStore = useChatStore();
const isUploading = ref(false);
let formData = new FormData();
const fileExtensionRegex = /\.(fasta|fastq|fna|fa|csv|tsv)$/i;

const canSelect = computed(()=>{
  const lastSearchMsgId = chatStore.currentChatMessages.filter(item=>item.agentType == AgentType.SEQUENCE_SEARCH).pop()?.id || '';
   return props.id == lastSearchMsgId && !isUploading.value
})

const canUpload = computed(()=>{
  let validate = true;
    selectedFiles.value.forEach((file:any)=> {
            if(!fileExtensionRegex.test(file.name)){
              validate = false;
            }
          });
   return validate && canSelect && !isUploading.value && !chatStore.currentChat.isGenerating
})

const selectedFiles = ref([] as File[]);
const rules = ref([
        (files:any) => {
          let validate = true;
          files.forEach((file:any)=> {
            if(!fileExtensionRegex.test(file.name)){
              validate = false;
            }
          });
          return validate||'You can only upload fasta/fastq/fna/fa/csv/tsv files';
        },
      ],);


async function handleUploadFiles() {
  isUploading.value = true;
  try {
    formData = new FormData();
    selectedFiles.value.forEach(file=>{
      formData.append('files', file);
    })
    const { data } = await uploadSequenceFiles(formData);
    if (data?.success == 'True') {
      // await chatStore.updateUploadSequenceFileMessage('', data.filenames, props.id)
      selectedFiles.value = [];
      formData = new FormData();
    }
  } catch (error) {
    console.error('upload file fail', error);
  } finally {
      isUploading.value = false;
  }
};



</script>

<style lang="scss" scoped>
.upload-sequence-file {
  min-width: 400px;
  .upload-section {
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

    .upload-btn {
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

    .loading-animation {

      // 简单的文字加载动画
      &::after {
        content: '.';
        animation: dots 1s steps(5, end) infinite;
      }

      // 定义动画效果
      @keyframes dots {

        0%,
        20% {
          color: rgba(0, 0, 0, 0);
          text-shadow:
            .25em 0 0 rgba(0, 0, 0, 0),
            .5em 0 0 rgba(0, 0, 0, 0);
        }

        40% {
          color: black;
          text-shadow:
            .25em 0 0 rgba(0, 0, 0, 0),
            .5em 0 0 rgba(0, 0, 0, 0);
        }

        60% {
          text-shadow:
            .25em 0 0 black,
            .5em 0 0 rgba(0, 0, 0, 0);
        }

        80%,
        100% {
          text-shadow:
            .25em 0 0 black,
            .5em 0 0 black;
        }
      }
    }
  }

  .upload-files-list {
    list-style: none;
    padding: 0;
    margin: 0;
    text-align: left;
    display: inline-block;
    width: 60%;
    margin: 10px;

    .upload-file {
      display: flex;
      justify-content: space-between;
      align-items: center;
      background-color: #444;
      padding: 5px 10px;
      border-radius: 5px;
      margin-bottom: 5px;
      box-shadow: 0 2px 4px rgba(0, 0, 0, 0.2);
      
    }
    
  }

}
</style>