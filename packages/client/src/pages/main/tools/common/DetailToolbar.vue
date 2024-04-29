<template>
  <div class='detail-toolbar'>
    <v-btn icon='mdi-arrow-left'
           :style="{color: '#ECECEC'}"
           variant='text'
           @click="$emit('onBack')"/>
    <div class='title'>
      {{title}}
    </div>
    <v-btn icon='mdi-dots-vertical'
           :style="{color: '#ECECEC'}"
           variant='text'/>
  </div>
</template>

<script setup lang='ts'>
const props = defineProps({
  title: String,
})
defineEmits([
  'onBack'
])
</script>

<style scoped lang='scss'>
.detail-toolbar {
  display: flex;
  align-items: center;
  justify-content: space-between;
  .title {
    color: #ECECEC;
    background-color: #1E2032;
    border-radius: 6px;
    font-size: 16px;
    height: 40px;
    line-height: 40px;
    padding: 0 15px;
  }
}
</style>