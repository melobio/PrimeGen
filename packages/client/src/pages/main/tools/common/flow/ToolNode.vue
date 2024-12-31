<template>
  <div class='tool-node'
    :style="borderStyle"
  >
    <div class='title'>
      {{localData.name}}
    </div>
    <div class='desc'>
      <span>
        {{localData.shortDesc || localData.desc}}
      </span>
    </div>
    <!--<NodeEditField-->
    <!--  required-->
    <!--  editable-->
    <!--  :label="'API Base'"-->
    <!--  class='mt-3'-->
    <!--  v-model='localData.apiBase'-->
    <!--/>-->
    <template
      v-for='(field) in localData.fields'
      :key='field.label'
    >
      <NodeEditField
        required
        editable
        :hide-edit-btn="!Boolean(localData.uuid)"
        :label="field.label"
        :placeholder="field.hint"
        class='mt-3'
        v-model='field.value'
      />
    </template>
    <NodeChildField
      v-for="(child, index) in localData.children"
      :ref="(el) => setChildrenRef(el, index)"
      :key="child.name"
      class='pl-5 mt-3'
      v-bind='child'
    />
    <!--<NodeOnlineField-->
    <!--  :enable-refresh="enableRefresh"-->
    <!--  @on-refresh="refreshState"-->
    <!--  :label="'Online State'"-->
    <!--  class='pl-5 pr-2 mt-3'-->
    <!--  v-model='localData.online'/>-->
    <div class='flex-1-1'/>
    <div class='title-bottom'
         :ref="(el) => setTargetRef(el)"
    >
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

  <Handle
    v-if="targetHandle"
    :key="targetHandle.uuid"
    :id="`t:${targetHandle.uuid}`"
    type="target"
    :position="Position.Right"
    :style="targetHandleStyle"
  />
</template>

<script setup lang='ts'>
import {computed, reactive, ref, watch} from 'vue'
import NodeEditField from '@/pages/main/tools/common/flow/fields/NodeEditField.vue'
import type {DeviceProp, ToolProp} from '@/pages/main/tools/common/flow/NodeTypes'
import {checkDeviceState, updateDevice, updateTool} from "@/api";
import NodeOnlineField from "@/pages/main/tools/common/flow/fields/NodeOnlineField.vue";
import {Handle, Position} from "@vue-flow/core";
import type {ToolNode} from "@xpcr/common";
import NodeChildField from "@/pages/main/tools/common/flow/fields/NodeChildField.vue";
import {NodeTypeColors} from "@xpcr/common/src/node";
import {useExperimentsStore} from "@/stores/experiments";

const props = defineProps<ToolProp>()
const expStore = useExperimentsStore();
const active = computed(() => expStore.activeNodeId === props.data.uuid);

const localData = reactive<ToolNode>(props.data);
const handles = ref(props.data.children?.map((child) => ({
  ...child,
  offsetTop: 0,
})));
const targetHandle = ref({
  offsetTop: 0,
  uuid: props.data.uuid,
});
const borderStyle = computed(() => {
  return {
    'border-color': active.value ? 'hsl(215 20.2% 65.1%)' : NodeTypeColors[localData.type],
  }
})
const targetHandleStyle = computed(() => {
  return {
    background: '#151728',
    width: '20px',
    height: '20px',
    borderRadius: '10px',
    border: `3px solid ${NodeTypeColors[localData.type]}`,
    boxSizing: 'border-box',
    right: '-8px',
    top: `${targetHandle.value.offsetTop}px`,
    bottom: 'auto',
  }
})

watch(() => localData, () => {
  // console.log('localData', localData);
  updateTool(localData.uuid, localData).then((res) => {
    if (res.success) {
      // update success
    }
  })
}, { deep: true })

function setChildrenRef(child: any, index: number) {
  if (!child) return;
  const el = child.$el as HTMLDivElement;
  if (handles.value) {
    handles.value[index].offsetTop = el.offsetTop + el.clientHeight / 2 + 3; // border width 3
  }
}
function setTargetRef(el: any) {
  if (!el) return;
  targetHandle.value.offsetTop = el.offsetTop + el.clientHeight / 2 + 3; // border width 3
}
</script>

<style scoped lang='scss'>
.tool-node {
  width: 348px;
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
    background-color: #151728;
    padding: 10px 15px;
    > span {
      line-height: 20px;
      text-align: start;
      color: #6D6F84;
      font-size: 14px;
      box-sizing: border-box;
      white-space: pre-wrap;
      display: -webkit-box;
      -webkit-box-orient: vertical;
      -webkit-line-clamp: 8;
      overflow: hidden;
    }
  }
  .title-bottom {
    color: #6D6F84;
    text-align: right;
    padding: 0 25px;
    margin-top: 20px;
    margin-bottom: 13px;
  }
}
</style>