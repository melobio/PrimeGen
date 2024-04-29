<template>
  <div class='experiment-container'>
    <ExpList />
    <ExpDetail
      v-if="experimentsStore.currentExperiment"
      v-bind="experimentsStore.currentExperiment"
      class='flex-1-1 ml-4'
    />
  </div>
</template>

<script setup lang='ts'>
import ExpList from '@/pages/main/experiment/ExpList.vue'
import ExpDetail from '@/pages/main/experiment/ExpDetail.vue'
import { onMounted } from 'vue'
import { useExperimentsStore } from '@/stores/experiments'
const experimentsStore = useExperimentsStore();

onMounted(() => {
  experimentsStore.getAllExperiments()
})
</script>

<style scoped lang='scss'>
.experiment-container {
  display: flex;
  flex-direction: row;
  align-items: stretch;
  width: 100%;
  height: 100%;
}
</style>