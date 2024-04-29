<template>
  <div class='tool-list-item' :style="{'box-shadow': `2px 0 0 0 ${item?.color}`}">
    <TintImg
      :src='item?.icon'
      :tint='color'
    />
    <span class='item-title' :style="{color: color}">
      {{item.name}}
    </span>
  </div>
</template>

<script setup lang='ts'>
import { ToolItem } from '@/pages/main/tools/ToolItem'
import TintImg from '@/components/TintImg.vue'
import { computed } from 'vue'

const props = defineProps({
  item: {
    type: ToolItem,
    default: () => ({}),
  },
  activeItem: {
    type: ToolItem,
    default: () => ({}),
  }
})

const color = computed(() => {
  if (props.activeItem.name === props.item.name) {
    return '#ECECEC'
  } else {
    return  '#6D6F84'
  }
})
</script>

<style scoped lang='scss'>
.tool-list-item {
  height: 48px;
  min-height: 48px;
  border-radius: 6px;
  display: flex;
  align-items: center;
  background-color: #1E2032;
  padding: 0 15px;
  cursor: pointer;
  .item-title {
    margin-left: 15px;
  }
}
.tool-list-item + .tool-list-item {
  margin-top: 5px;
}
</style>