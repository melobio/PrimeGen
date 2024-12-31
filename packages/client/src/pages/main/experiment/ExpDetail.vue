<template>
  <div class='exp-detail'>
    <VueFlow
      v-if="elements"
      ref='flow'
      v-model="elements"
      :node-types='NodeTypes'
      class='w-100 h-100'
      fit-view-on-init
      :max-zoom="1"
      :default-viewport="{ zoom: 1 }"
    >
      <Background :pattern-color="'#6D6F84'" :gap="20" />
    </VueFlow>
  </div>
</template>

<script setup lang='ts'>
import '@vue-flow/core/dist/style.css';
import '@vue-flow/core/dist/theme-default.css';
import '@vue-flow/controls/dist/style.css';
import '@vue-flow/minimap/dist/style.css';
import '@vue-flow/node-resizer/dist/style.css'
import {VueFlow, useVueFlow } from '@vue-flow/core'
import type { Elements } from '@vue-flow/core'
import { Background } from '@vue-flow/background'
import { NodeTypes } from '@/pages/main/tools/common/flow/NodeTypes'
import {nextTick, onMounted, ref, watch} from "vue";
import {commonToolInfo} from "@/api";
import type {ExperimentInterface, NodeChild} from "@xpcr/common";
import { NodeTypes as NT } from '@xpcr/common/src/node';
import type { Node } from '@xpcr/common/src/node';
import {useExperimentsStore} from "@/stores/experiments";
const props = defineProps<ExperimentInterface>()
const { onNodeClick } = useVueFlow();
const expStore = useExperimentsStore();

const flow = ref()
const elements = ref<Elements>();

let position: {x: number, y: number} = { x: 0, y: 0 };

function getNextPosition() {
  // position.x -= 400;
  const ret = {...position};
  position.y += 500;
  if (position.y > 1000) {
    position.x -= 500;
    position.y = 0;
  }
  return ret;
}

async function loadChildren(els: Elements, parentUUID: string, children: NodeChild[]) {
  const grandson: { [key: string]: NodeChild[] } = {};
  for (const { type, uuid } of children) {
    if (uuid && type) {
      const { success, data } = await commonToolInfo<Node>(uuid, type);
      if (success) {
        // node
        const nodeExists = els.findIndex((el) => el.id === uuid) > -1;
        if (!nodeExists) {
          els.push({
            id: uuid,
            type: type,
            position: getNextPosition(),
            data: data,
          })
        }
        // link
        els.push({
          id: `${parentUUID}:${uuid}`,
          source: parentUUID,
          target: uuid,
          animated: false,
          sourceHandle: `s:${uuid}`,
          targetHandle: `t:${uuid}`,
          style: () => ({
            stroke: '#1484FC',
            strokeWidth: 2,
          }),
        })
        // sub node
        grandson[uuid] = data?.children || [];
      }
    }
  }

  for (const uuid in grandson) {
    const children = grandson[uuid];
    await loadChildren(els, uuid, children);
  }
}

onMounted(async () => {
  const els: Elements = [
    {
      id: props.uuid,
      type: NT.Experiment,
      position: { ...position },
      data: props,
    },
  ];
  position.x -= 500;
  await loadChildren(els, props.uuid, props.nodeChildren || [])
  elements.value = els;
});
watch(() => props.uuid, async () => {
  elements.value = [];
  await nextTick(async () => {
    position.x = 0;
    position.y = 0;
    const els: Elements = [
      {
        id: props.uuid,
        type: NT.Experiment,
        position: { ...position },
        data: props,
      },
    ];
    position.x -= 500;
    await loadChildren(els, props.uuid, props.nodeChildren || [])
    elements.value = els;
  });
}, {deep: true});

onNodeClick((ev) => {
  // console.log('click', ev);
  expStore.setActiveNodeId(ev.node.id);
})
</script>

<style scoped lang='scss'>
.exp-detail {
  height: 100%;
  width: 100%;
}
</style>