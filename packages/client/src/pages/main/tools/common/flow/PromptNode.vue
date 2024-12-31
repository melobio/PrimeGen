<template>
  <div class='prompt-node'
    :style="borderStyle"
  >
    <div class='title'>
      {{localData.name}}
    </div>
    <div class='desc'>
      {{localData.shortDesc}}
    </div>
    <template
      v-for='(field) in localData.fields'
      :key='field.label'
    >
      <NodeEditField
        v-if="field.type === NodeFieldType.TextField"
        required
        editable
        :hide-edit-btn="!Boolean(localData.uuid)"
        :label="field.label"
        :placeholder="field.hint"
        class='mt-3'
        v-model='field.value'
      />
      <NodeTextAreaField
        v-else-if="field.type === NodeFieldType.TextArea"
        required
        editable
        :hide-edit-btn="!Boolean(localData.uuid)"
        :label="field.label"
        :placeholder="field.hint"
        class='mt-3'
        v-model='field.value'
      />
    </template>
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

    <!--<PromptTempEditDialog-->
    <!--  v-model='templateEditDialog'-->
    <!--  v-model:template='localData.template'-->
    <!--  :variable='variable'-->
    <!--/>-->
  </div>

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
import { computed, isReactive, onMounted, reactive, ref, watch } from 'vue'
import NodeEditField from '@/pages/main/tools/common/flow/fields/NodeEditField.vue'
import { updatePrompt } from '@/api'
import type { PromptProp } from '@/pages/main/tools/common/flow/NodeTypes'
import type {PromptNode} from "@xpcr/common";
import {Handle, Position} from "@vue-flow/core";
import {NodeFieldType} from "@xpcr/common/src/node-field";
import NodeTextAreaField from "@/pages/main/tools/common/flow/fields/NodeTextAreaField.vue";
import {NodeTypeColors} from "@xpcr/common/src/node";
import {useExperimentsStore} from "@/stores/experiments";
const props = defineProps<PromptProp>()

const localData = reactive<PromptNode>(props.data);

const expStore = useExperimentsStore();
const active = computed(() => expStore.activeNodeId === props.data.uuid);

watch(() => localData, () => {
  // console.log('localData', localData);
  if (localData.uuid) {
    updatePrompt(localData.uuid, localData).then((res) => {
      if (res.success) {
        // update success
      }
    })
  }
}, { deep: true })

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
function setTargetRef(el: any) {
  if (!el) return;
  targetHandle.value.offsetTop = el.offsetTop + el.clientHeight / 2 + 3; // border width 3
}
</script>

<style scoped lang='scss'>
.prompt-node {
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
    line-height: 40px;
    height: 40px;
    text-align: start;
    padding: 0 15px;
    color: #6D6F84;
    font-size: 14px;
    overflow: hidden;
    text-overflow: ellipsis;
    white-space: nowrap;
    background-color: #151728;
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