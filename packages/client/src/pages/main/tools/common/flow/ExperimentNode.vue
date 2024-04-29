<template>
  <div class='experiment-node'
    :style="borderStyle"
  >
    <div class='title'>
      {{localData.name}}
    </div>
    <div class='desc'>
      {{localData.desc}}
    </div>
    <!--<NodeEditField-->
    <!--  required-->
    <!--  editable-->
    <!--  :label="'API Base'"-->
    <!--  class='pl-5 pr-2 mt-3'-->
    <!--  v-model='localData.apiBase'-->
    <!--/>-->
    <!--<NodeFuncField-->
    <!--  v-for="(func) in localData.functions"-->
    <!--  :key="func.name"-->
    <!--  class='pl-5 mt-3'-->
    <!--  v-bind='func'-->
    <!--/>-->
    <NodeChildField
      v-for="(child, index) in localData.nodeChildren"
      :ref="(el) => setChildrenRef(el, index)"
      :key="child.name"
      class='pl-5 mt-3'
      v-bind='child'
    />
    <div class='flex-1-1'/>
    <div class='title-bottom'>
      {{localData.name}}
    </div>
    <div class='button-bar' v-if="active">
      <v-btn
        width='40px'
        height='40px'
        icon='mdi-trash-can-outline'
        color='#292B3C'
        rounded
        :style="{color: '#6D6F84'}"
        :elevation='0'/>
      <v-btn
        width='40px'
        height='40px'
        icon='mdi-content-copy'
        color='#292B3C'
        rounded
        class='ml-1'
        :style="{color: '#6D6F84'}"
        :elevation='0'/>
      <v-btn
        width='40px'
        height='40px'
        icon='mdi-swap-horizontal'
        color='#292B3C'
        rounded
        class='ml-1'
        :style="{color: '#6D6F84'}"
        :elevation='0'/>
    </div>

  </div>

  <Handle
    v-for="(child) in handles"
    :key="child.uuid"
    :id="`s:${child.uuid}`"
    type="source"
    :position="Position.Left"
    :style="{
      background: '#151728',
      width: '20px',
      height: '20px',
      borderRadius: '10px',
      border: '2px solid #1484FC',
      boxSizing: 'border-box',
      left: '-8px',
      top: `${child.offsetTop}px`,
      bottom: 'auto',
    }"
  />
</template>

<script setup lang='ts'>
import { Handle, Position } from '@vue-flow/core'
import {computed, ref} from 'vue'
import type {ExperimentProp} from '@/pages/main/tools/common/flow/NodeTypes'
import NodeChildField from "@/pages/main/tools/common/flow/fields/NodeChildField.vue";
import {NodeTypeColors} from "@xpcr/common/src/node";
import type {ExperimentInterface} from "@xpcr/common";
import { NodeTypes } from '@xpcr/common/src/node';
import {useExperimentsStore} from "@/stores/experiments";

const props = defineProps<ExperimentProp>()

const expStore = useExperimentsStore();
const active = computed(() => expStore.activeNodeId === props.data.uuid);

const localData = ref<ExperimentInterface>(props.data);
const handles = ref(props.data.nodeChildren?.map((child) => ({
  ...child,
  offsetTop: 0,
})));

const borderStyle = computed(() => {
  return {
    'border-color': active.value ? 'hsl(215 20.2% 65.1%)' : NodeTypeColors[NodeTypes.Experiment],
  }
})
function setChildrenRef(child: any, index: number) {
  if (!child) return;
  const el = child.$el as HTMLDivElement;
  if (handles.value) {
    handles.value[index].offsetTop = el.offsetTop + el.clientHeight / 2 + 3; // border width 3
  }
}
</script>

<style scoped lang='scss'>
.experiment-node {
  width: 348px;
  //height: 382px;
  background-color: #1E2032;
  border: 3px solid #1484FC;
  border-radius: 24px;
  position: relative;
  display: flex;
  flex-direction: column;
  .button-bar {
    position: absolute;
    top: -48px;
    left: 0;
    right: 0;
    display: flex;
    align-items: center;
    justify-content: center;
  }
  .title {
    font-size: 18px;
    color: #cecece;
    line-height: 45px;
    height: 45px;
    text-align: center;
  }
  .desc {
    line-height: 20px;
    //height: 40px;
    text-align: start;
    padding: 10px 15px;
    color: #6D6F84;
    font-size: 14px;
    //overflow: hidden;
    //text-overflow: ellipsis;
    //white-space: nowrap;
    background-color: #151728;
    white-space: pre-wrap;
  }
  .title-bottom {
    color: #6D6F84;
    text-align: right;
    padding: 0 25px;
    line-height: 25px;
    height: 25px;
    margin-top: 20px;
    margin-bottom: 15px;
  }
}
</style>