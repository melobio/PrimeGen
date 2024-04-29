<template>
  <div class='tools-container'>
    <ToolList />
    <!--<Agents v-if="toolsStore.currentTool.name === 'Agents'"/>-->
    <!--<Inputs v-if="toolsStore.currentTool.name === 'Input'"/>-->
    <!--<Prompts v-if="toolsStore.currentTool.name === 'Prompts'"/>-->
    <CardList />
  </div>
</template>

<script setup lang='ts'>
import ToolList from '@/pages/main/tools/ToolList.vue'

import { useToolsStore } from '@/stores/tools'
import CardList from '@/pages/main/tools/CardList.vue'
const toolsStore = useToolsStore();
</script>

<style scoped lang='scss'>
.tools-container {
  flex: 1;
  display: flex;
  flex-direction: row;
  align-items: stretch;
  width: 100%;
  height: 100%;
}
</style>