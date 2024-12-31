<template>
    <v3-drag-zoom-container
    >
        <div class="player" ref="dplayer" style="width:400px;height:200px"  ></div>

  </v3-drag-zoom-container>
</template>

<script setup lang="ts">

// import DPlayer from 'dplayer';
import { onMounted, ref } from 'vue';
// import video from '@/assets/demo2.mp4';
import {V3DragZoomContainer} from "v3-drag-zoom";

const dplayer = ref(null)

// onMounted(() => {
//   new DPlayer({
//     container: dplayer.value,
//     video: {
//       url: video,
//       controller: true,
//     }
//   })
// })
const dragRef = ref();

const onChange = () => {
  setTimeout(() => {
    dragRef.value.reset();
  });
};


</script>