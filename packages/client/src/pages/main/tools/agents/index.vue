<!--<template>-->
<!--  <div class='agents-container'>-->
<!--    <div class='flex-1-1 d-flex flex-column overflow-hidden'>-->
<!--      <div class='text-center'>-->
<!--        <div class='tools-sub-count mt-10 mb-5 d-inline-block'>-->
<!--          Agents({{toolsStore.agents.length}})-->
<!--        </div>-->
<!--      </div>-->
<!--      <div class='flex-1-1 h-100 d-flex flex-row flex-wrap overflow-y-auto align-content-start'>-->
<!--        <AgentItem-->
<!--          v-for='(item) in toolsStore.agents'-->
<!--          :key='item.id'-->
<!--          :item='item'-->
<!--        />-->
<!--      </div>-->
<!--    </div>-->
<!--  </div>-->
<!--</template>-->

<!--<script setup lang='ts'>-->
<!--import { useToolsStore } from '@/stores/tools'-->
<!--import AgentItem from '@/pages/main/tools/agents/AgentItem.vue'-->
<!--import { onMounted } from 'vue'-->
<!--import { toolsAgents } from '@/api'-->
<!--import type { Agent } from '@/stores/tools-types'-->
<!--const toolsStore = useToolsStore();-->
<!--onMounted(async () => {-->
<!--  const { success, data } =  await toolsAgents<Agent[]>()-->
<!--  if (success) {-->
<!--    toolsStore.agents = data;-->
<!--  }-->
<!--})-->
<!--</script>-->

<!--<style scoped lang='scss'>-->
<!--.agents-container {-->
<!--  width: 100%;-->
<!--  height: 100%;-->
<!--  flex: 1;-->
<!--  display: flex;-->
<!--  flex-direction: row;-->
<!--  align-content: start;-->
<!--  .tools-sub-count {-->
<!--    height: 40px;-->
<!--    width: 240px;-->
<!--    background-color: #1E2032;-->
<!--    color: #ECECEC;-->
<!--    border-radius: 6px;-->
<!--    line-height: 40px;-->
<!--    text-align: center;-->
<!--  }-->
<!--}-->
<!--</style>-->