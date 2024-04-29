<template>
  <div class='exp-list'>
    <div class='search-container'>
      <div class='search'>
        <v-icon icon='mdi-magnify' style='color: #545665'/>
        <input class='search-input' placeholder='Search for experiment'/>
      </div>
      <v-btn icon="mdi-plus" size='40px' class='ml-3' rounded color='#1E2032'></v-btn>
    </div>
    <div class='list-container'>
      <exp-list-item
        v-for='(exp) in experimentsStore.experiments'
        v-bind='exp'
        :key='exp.uuid'
        @click="select(exp)"
      />
    </div>
  </div>
</template>

<script setup lang='ts'>
import ExpListItem from '@/pages/main/experiment/ExpListItem.vue'
import { useExperimentsStore } from '@/stores/experiments'
import type {ExperimentInterface} from "@xpcr/common/src/experiment";
const experimentsStore = useExperimentsStore();
function select(exp: ExperimentInterface) {
  experimentsStore.setCurrentExperiment(exp);
}
</script>

<style scoped lang='scss'>
.exp-list {
  width: 300px;
  min-width: 300px;
  height: 100%;
  padding: 15px 0;
  display: flex;
  flex-direction: column;
  .search-container {
    //padding: 0 15px;
    padding-left: 15px;
    margin-top: 25px;
    display: flex;
    flex-direction: row;
    .search {
      flex: 1;
      height: 40px;
      background: #1E2032;
      border-radius: 6px;
      display: flex;
      align-items: center;
      padding: 0 15px;
      .search-input {
        height: 100%;
        flex: 1;
        margin-left: 10px;
        border: none;
        outline: none;
        color: white;
        &::placeholder {
          color: #545665;
        }
      }
    }
  }

  .list-container {
    padding: 20px 0 0 15px;
    overflow-y: auto;
    flex: 1;
  }
}
</style>