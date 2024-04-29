<template>
  <div class='img-container'
    :style="{
      width: `${width}px`,
      height: `${height}px`,
    }">
    <img
      :src='src' alt=''
      :style='{
         filter: `drop-shadow(${tint} 40px 0)`
      }'
    />
  </div>
</template>

<script setup lang='ts'>
const props = defineProps({
  width:{
    type: Number,
    default: 20,
  },
  height: {
    type: Number,
    default: 20,
  },
  src: {
    type: String,
  },
  tint: {
    type: String,
    default: '#31DA9F',
  }
})
</script>

<style scoped lang='scss'>
.img-container {
  overflow: hidden;
  //position: relative;
  display: inline-block;
  img {
    width: 100%;
    height: 100%;
    margin-left: -40px;
  }
}
</style>