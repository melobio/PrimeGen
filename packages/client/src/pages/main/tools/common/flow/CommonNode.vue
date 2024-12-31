<template>
  <div class='common-node'>
    <div class='title'>
      {{localData.name}}
    </div>
    <div class='desc'>
      {{localData.shortDesc || localData.desc}}
    </div>
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
</template>

<script setup lang='ts'>
import {computed, ref, watch} from 'vue'
import type { CommonProp } from '@/pages/main/tools/common/flow/NodeTypes'
import { updateCommonToolItem } from '@/api'
import type { Node } from '@xpcr/common'
import {useExperimentsStore} from "@/stores/experiments";
const expStore = useExperimentsStore();

const props = defineProps<CommonProp>()

const localData = ref<Node>(props.data);

watch(() => localData, () => {
  // console.log('localData', localData);
  updateCommonToolItem(localData.value.uuid, props.type, localData.value).then((res) => {
    if (res.success) {
      // update success
    }
  })
}, { deep: true })
const active = computed(() => expStore.activeNodeId === props.data.uuid);
</script>

<style scoped lang='scss'>
.common-node {
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
    margin-top: 20px;
    margin-bottom: 13px;
  }
}
</style>