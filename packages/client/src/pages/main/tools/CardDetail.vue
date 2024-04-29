<template>
  <div class='card-detail'>
    <VueFlow
      v-if="elements"
      ref='flow'
      v-model="elements"
      :node-types='NodeTypes'
      class='position-absolute w-100 h-100'
      fit-view-on-init
      :max-zoom="1"
      :default-viewport="{ zoom: 1 }"
    >
      <Background :pattern-color="'#6D6F84'" :gap="20" />
    </VueFlow>

    <DetailToolbar
      class='position-absolute w-100 mt-10'
      :title="`${item?.name} Detail Information`"
      @onBack="$emit('onBack')"
    />
  </div>
</template>

<script setup lang='ts'>
import '@vue-flow/core/dist/style.css';
import '@vue-flow/core/dist/theme-default.css';
import '@vue-flow/controls/dist/style.css';
import '@vue-flow/minimap/dist/style.css';
import '@vue-flow/node-resizer/dist/style.css'
import { Panel, VueFlow, isNode, useVueFlow } from '@vue-flow/core'
import DetailToolbar from '@/pages/main/tools/common/DetailToolbar.vue'
import { Background } from '@vue-flow/background'
import {onMounted, ref} from 'vue'
import type {PropType} from 'vue'
import { NodeTypes } from '@/pages/main/tools/common/flow/NodeTypes'
import type { Node } from '@xpcr/common/src/node'
import { useToolsStore } from '@/stores/tools'
import type { Elements } from '@vue-flow/core'
import {commonToolInfo} from "@/api";

const toolsStore = useToolsStore();

const props = defineProps({
  item: {
    type:  Object as PropType<Node>,
    default: () =>({}),
  }
})

defineEmits([
  'onBack'
])
const flow = ref()
const elements = ref<Elements>();

async function loadChildren() {
  const eles: Elements = [];
  for (const { type, uuid } of (props.item.children || [])) {
    if (uuid && type) {
      const { success, data } = await commonToolInfo(uuid, type);
      if (success) {
        // console.log('data',data);
        eles.push({
          id: uuid,
          type: type,
          position: { x: -600, y: 100 },
          data: data,
        })
        eles.push({
          id: `${props.item.uuid}:${uuid}`,
          source: props.item.uuid,
          target: uuid,
          animated: true,
          sourceHandle: `s:${uuid}`,
          targetHandle: `t:${uuid}`,
          style: () => ({
            stroke: '#1484FC',
            strokeWidth: 2,
          }),
        })
      }
    }
  }
  return eles;
}

onMounted(async () => {
  const eles: Elements = [
    {
      id: props.item.name,
      type: toolsStore.currentTool.name,
      position: { x: 0, y: 0 },
      data: props.item,
    },
  ];
  eles.push(...await loadChildren());
  // console.log('eles', eles);
  elements.value = eles;
});

</script>

<style scoped lang='scss'>
.card-detail {
  width: 100%;
  height: 100%;
  flex: 1;
  position: relative;
}
</style>