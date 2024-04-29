<script setup lang="ts">
import {computed} from "vue";

const props = defineProps({
  data: String,
  plot: {
    type: Array,
    default: () => ([]),
  },
  success: Boolean,
  imgWidth: {
    type: Number,
    default: 500,
  }
})

const picSize = {
  width: 1280,
  height: 720,
}

const thePlot: any[] = props.plot.length > 0 ? props.plot : [
  0,
  0,
  0,
  0,
  0,
  0,
]

const boundingBoxStyle = computed(() => {
  const scale = props.imgWidth / picSize.width
  return {
    left: `${thePlot[0] * scale}px`,
    top: `${thePlot[1] * scale}px`,
    width: `${(thePlot[2] - thePlot[0]) * scale}px`,
    height: `${(thePlot[3] - thePlot[1]) * scale}px`,
    display: thePlot[4] ? 'block' : 'none',
  }
})

const tipTagStyle = computed(() => {
  const scale = props.imgWidth / picSize.width
  return {
    left: `${thePlot[0] * scale}px`,
    top: `${thePlot[1] * scale - 26}px`,
    display: thePlot[4] ? 'block' : 'none',
  }
})

const detectTip = computed(() => {
  return thePlot[4] !== 0
})

</script>

<template>
  <div class="jetson-fault-check"
       :style="{ width: `${imgWidth}px`, 'border-color': success ? 'green' : 'red' }">
    <img :src="'data:image/png;base64,'+data" alt=""/>
    <div class="bounding-box" :style="boundingBoxStyle"/>
    <div class="tip-tag" :style="tipTagStyle">tip {{ parseFloat(thePlot[4]).toFixed(2) }}</div>
    <span
      style="padding-left: 10px; font-size: 0.8rem"
      :style="{ color: success ? 'green' : 'red'  }"
    >
      {{ detectTip ? 'Has tips' : 'No tips' }}
    </span>
  </div>
</template>

<style scoped lang="scss">
.jetson-fault-check {
  position: relative;
  border: 1px solid transparent;
  img {
    width: 100%;
    height: 100%;
    display: inline-block;
  }
  .bounding-box {
    position:absolute;
    //background-color: red;
    border: 1px solid red;
  }
  .tip-tag {
    white-space: pre;
    position:absolute;
    background-color: white;
    color: red;
    opacity: 0.8;
    font-weight: bold;
    height: 25px;
    line-height: 25px;
    display: flex;
    align-items: center;
    padding: 0 3px;
    font-size: 12px;
  }
}
</style>