<template>
  <div class="audio-record-button">
    <v-btn
      :disabled="chatStore.currentChat.isGenerating"
      elevation='0' size='60' color='#292B3C' @click="toggleRecording">
      <img v-if="!chatStore.isRecordingAudio" src='@/assets/audio-record.svg' alt='' />
      <img  v-else src='@/assets/audio-record-loading.svg' alt=''/>
    </v-btn>
  </div>
</template>
<script setup lang="ts">
import { onUnmounted, ref } from "vue";
import { useChatStore } from '@/stores/chats';
defineProps({
  disabled:{
    type:Boolean,
    default:false
  }
})

const mediaRecorder = ref<MediaRecorder | any>()
const chunks = ref<Array<Blob>>([])
const recordedAudioURL = ref('')
let stream: MediaStream
const chatStore = useChatStore();

onUnmounted(() => {
  closeStream(stream)
})

const initRecording = () => {
  navigator.mediaDevices.getUserMedia({ audio: true })
    .then(audioStream => {
      stream = audioStream;
      mediaRecorder.value = new MediaRecorder(audioStream as MediaStream, {
        mimeType: 'audio/webm',
      });
      mediaRecorder.value.ondataavailable = (e: BlobEvent) => {
        if (e.data && e.data.size > 0) {
          chunks.value.push(e.data);
        }
      };
      mediaRecorder.value.onstop = () => {
        const audioBlob = new Blob(chunks.value, { type: 'audio/webm' });
        recordedAudioURL.value = URL.createObjectURL(audioBlob);
        console.log('recordedAudioURL===>', recordedAudioURL)
        convertToBase64(audioBlob)
        stream.getTracks().forEach((track) => track.stop());
      };
      mediaRecorder.value.start();
      chatStore.setRecordingAudioState(true)
    })
    .catch(error => {
      console.error('record error:', error);
    });
}

const convertToBase64 = (audioBlob: Blob) => {
  // if (chunks.value.length > 0) {
  //   const reader = new FileReader();
  //   reader.onloadend = () => {
  //     const base64DataString = reader.result as string;
  //     console.log('Audio===>base64 data:', base64DataString);
  //   };
  //   reader.readAsDataURL(audioBlob);
  // }
  chatStore.sendVoice(audioBlob)
}
const toggleRecording = () => {
  if (chatStore.isRecordingAudio) {
    mediaRecorder.value.stop()
    chatStore.setRecordingAudioState(false)
    chunks.value = []
  } else {
    initRecording()
  }
}

const closeStream = (stream: MediaStream) => {
  if (stream) {
    console.log(stream.getTracks())
    stream.getTracks().forEach((track) => track.stop())
  }
}



</script>
<style scoped lang="scss">
.audio-record-button {
  button {
    border-radius: 10px;
    margin-right: 10px;
  }

  img {
    width: 25px;
    height: 25px;
  }

  svg path,
  svg rect {
    fill: #31DA9F;
  }
}
</style>