<template>
  <div class='agent-item'>
    <div class='agent-title'>
      {{item.name}}
    </div>
    <div class='agent-desc'>
      {{item.desc}}
    </div>
    <div class='agent-version'>
      <div class='version'>
        v{{item.version}}
      </div>
      <v-btn variant='text'>
        <span class='btn-detail'>
          Detail View
        </span>
      </v-btn>
    </div>
  </div>
</template>

<script setup lang='ts'>
import type { AgentNode } from '@xpcr/common/src/agent'
import type { PropType } from 'vue'

defineProps({
  item: {
    type: Object as PropType<AgentNode>,
    default: () =>({}),
  }
})
</script>

<style scoped lang='scss'>

.agent-item {
  width: 290px;
  //height: 162px;
  background-color: #292B3C;
  margin-bottom: 15px;
  margin-right: 15px;
  border-radius: 10px;
  display: flex;
  flex-direction: column;
  .agent-title {
    color: #ECECEC;
    font-size: 16px;
    margin: 15px 15px 5px 15px;
  }
  .agent-desc {
    color: #545665;
    margin-bottom: 15px;
    margin-left: 15px;
    margin-right: 15px;
    flex: 1;
  }
  .agent-version {
    border-top: 1px solid #151728;
    display: flex;
    height: 40px;
    align-items: center;
    justify-content: space-between;
    padding: 0 0 0 15px;
    .version {
      color: #6D6F84;
    }
    .btn-detail {
      color: #31DA9F;
      text-transform: none;
    }
  }
}
</style>