<template>
  <div class='custom-input' :class="{active}">
    <slot name='icon'></slot>
    <input
      @keydown.enter='onEnter'
      :value='modelValue'
      @input='input'
      :placeholder='placeholder'
      @focus='active = true' @blur='active = false' :type="password ? 'password' : ''">
  </div>
</template>

<script setup lang='ts'>
import { ref } from 'vue'
// import { VIcon } from "vuetify"

const props = defineProps({
  modelValue: {
    type: String,
  },
  password: {
    type: Boolean,
    default: false,
  },
  placeholder: String,
})
const emit = defineEmits(['update:modelValue', 'onEnter'])

const active = ref(false)

function input(ev: Event) {
  const value = (ev.target as HTMLInputElement).value;
  // console.log(value);
  emit('update:modelValue', value);
}
function onEnter() {
  // console.log('onEnter');
  emit('onEnter')
}
</script>

<style scoped lang='scss'>
.custom-input {
  background: #151728;
  border-radius: 5px;
  height: 48px;
  min-height: 48px;
  border: 1px solid #151728;
  display: flex;
  flex-direction: row;
  align-items: center;
  padding: 0 10px;
  &:hover, &.active {
    border-color: #7b67ee;
  }
  input {
    width: 100%;
    height: 100%;
    border: none;
    outline: none;
    color: #ECECEC;
    margin-left: 10px;
  }
}
</style>