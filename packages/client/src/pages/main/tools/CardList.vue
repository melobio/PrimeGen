<template>
  <div class='w-100 h-100 flex-1-1 d-flex align-content-start'>
    <!--Card List-->
    <div class='flex-1-1 d-flex flex-column overflow-hidden' v-if='!isDetail'>
      <!--Tools Type and Count-->
      <div class='text-center'>
        <div class='tools-sub-count mt-10 mb-5 d-inline-block'>
          <!--Prompts({{toolsStore.prompts.length}})-->
          {{toolsStore.currentTool.name}}({{toolsStore.currentToolList.value.length}})
        </div>
      </div>
      <!--Card Flex List-->
      <div class='d-flex flex-row flex-wrap flex-1-1 overflow-y-auto align-content-start'>
        <CardItem
          v-for='(item) in toolsStore.currentToolList.value'
          :key='item.id'
          :item='item'
          @onViewDetail='viewDetail(item)'
        />
      </div>
    </div>
    <!--Card Detail-->
    <div class='flex-1-1 d-flex flex-row' v-else>
      <CardDetail
        :item='detailItem'
        @onBack="viewList"
      />
    </div>
  </div>
</template>

<script setup lang='ts'>
import { ref, watch } from 'vue'
import type { Node } from '@xpcr/common'
import { useToolsStore } from '@/stores/tools'
import CardItem from '@/pages/main/tools/CardItem.vue'
import CardDetail from '@/pages/main/tools/CardDetail.vue'
const toolsStore = useToolsStore();
const isDetail = ref(false);
const detailItem = ref<Node>()

watch(() => toolsStore.currentTool, async () => {
  viewList();
  await toolsStore.getCurrentToolList()
}, { immediate: true });

function viewDetail(item: Node) {
  isDetail.value = true;
  detailItem.value = item;
}
function viewList() {
  isDetail.value = false;
}
</script>

<style scoped lang='scss'>
.tools-sub-count {
  height: 40px;
  width: 240px;
  background-color: #1E2032;
  color: #ECECEC;
  border-radius: 6px;
  line-height: 40px;
  text-align: center;
}
</style>