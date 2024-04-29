<template>
  <div class="node-child-field">
    <span class="type">
      {{name}}
    </span>
    <v-icon icon='mdi-asterisk' size='12' color='red' v-if='required'/>
  </div>
</template>

<script setup lang='ts'>
import type {NodeChild} from "@xpcr/common/dist/node";

// defineProps({
//   type: String,
//   required: {
//     type: Boolean,
//     default: false,
//   },
// })
defineProps<NodeChild>()
// console.log(props);
</script>

<style scoped lang='scss'>
.node-child-field {
  height: 40px;
  min-height: 40px;
  flex-shrink: 0;
  align-items: center;
  box-sizing: border-box;
  display: flex;
  .type {
    text-transform: capitalize;
    color: #6D6F84;
    margin-right: 5px;
  }
}
</style>