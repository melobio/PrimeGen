<template>
  <div class='exp-list-item'>
    <div class='item-title'>
      {{name}}
    </div>
    <div class='item-content'>
      {{desc}}
    </div>
    <div class='item-by'>
      <div>
        by {{creatorName}}
      </div>
      <div>
        {{createTime}}
      </div>
    </div>
    <div class='item-states' v-if="expStore.currentExperiment?.uuid === uuid">
      <v-icon icon='mdi-check' color='#31DA9F'/>
    </div>
  </div>
</template>

<script setup lang='ts'>
import type { ExperimentInterface } from '@xpcr/common/src/experiment'
import {useExperimentsStore} from "@/stores/experiments";

const props = defineProps<ExperimentInterface>()
const expStore = useExperimentsStore()
</script>

<style scoped lang='scss'>
.exp-list-item {
  border-radius: 10px;
  background-color: #1E2032;
  min-height: 100px;
  padding: 15px;
  transition: background-color 0.3s linear;
  cursor: pointer;
  display: flex;
  flex-direction: column;
  overflow: hidden;
  position: relative;
  &:hover {
    background-color: #292B3C;
  }
  .item-title {
    color: #ECECEC;
    font-size: 16px;
    text-overflow: ellipsis;
    overflow: hidden;
    white-space: nowrap;
    height: 23px;
    line-height: 23px;
  }
  .item-content {
    margin-top: 5px;
    color: #545665;
    flex: 1;
    overflow: hidden;
    text-overflow: ellipsis;
    display: -webkit-box;
    -webkit-line-clamp: 4;
    -webkit-box-orient: vertical;
  }
  .item-by {
    margin-top: 13px;
    height: 20px;
    line-height: 20px;
    display: flex;
    align-items: center;
    color: #545665;
    justify-content: space-between;
  }
  .item-states {
    position: absolute;
    top: 10px;
    right: 15px;
  }
}
.exp-list-item + .exp-list-item {
  margin-top: 10px;
}
</style>