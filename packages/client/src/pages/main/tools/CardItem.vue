<template>
  <div class='card-item'>
    <div class='card-title'>
      {{item.name}}
    </div>
    <div class='card-desc'>
      <div class="desc">
        {{item.desc}}
      </div>
    </div>
    <div class='card-version'>
      <div class='version'>
        v{{item.version}}
      </div>
      <v-btn variant='text' @click="$emit('onViewDetail')">
        <span class='btn-detail'>
          Detail View
        </span>
      </v-btn>
    </div>
  </div>
</template>

<script setup lang='ts'>
import type { Node } from '@xpcr/common/src/node'
import type { PropType } from 'vue'

defineProps({
  item: {
    type: Object as PropType<Node>,
    default: () =>({}),
  }
})

defineEmits([
  "onViewDetail"
])
</script>

<style scoped lang='scss'>

.card-item {
  width: 290px;
  //height: 162px;
  background-color: #292B3C;
  margin-bottom: 15px;
  margin-right: 15px;
  border-radius: 10px;
  display: flex;
  flex-direction: column;
  .card-title {
    color: #ECECEC;
    font-size: 16px;
    margin: 15px 15px 5px 15px;
  }
  .card-desc {
    color: #545665;
    margin-bottom: 15px;
    margin-left: 15px;
    margin-right: 15px;
    line-height: 20px;
    flex: 1;
    overflow: hidden;
    display: flex;
    flex-direction: column;
    max-height: 80px;
    .desc {
      flex: 1;
    }
  }
  .card-version {
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