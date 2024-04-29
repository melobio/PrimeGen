<template>
  <div class='node-inputs-field'
    :style="{'padding-right': hideEditBtn ? '15px' : '0px'}"
  >
    <div class='label-container'>
      <span class='label'>
        {{label}}
      </span>
      <v-icon icon='mdi-asterisk' size='12' color='red' v-if='required'/>
    </div>
    <div
      class='value-container'
      v-for="(input, index) in localValue"
      :key="index"
    >
      <div
        :class="{active: active === index}"
        class='value'
        @mousedown.stop='() => {}'
      >
        <input
          :value="input"
          @change="(e) => localValue[index] = (e.target as HTMLInputElement).value"
          :readonly='readonly'
          :style="{'pointer-events': readonly ? 'none' : 'auto'}"
          :placeholder="placeholder"
          :type='inputType'
          @focus="active = index"
          @blur="active = -1"
        />
      </div>
      <template v-if="!hideEditBtn">
        <v-btn
          v-if="index === localValue.length - 1"
          variant='text'
          :style="{color: '#6D6F84'}"
          icon='mdi-plus'
          @click="addInputs"
        />
        <v-btn
          v-else
          variant='text'
          :style="{color: '#6D6F84'}"
          icon='mdi-close'
          @click="removeInputs(index)"
        />
      </template>
    </div>
  </div>
</template>

<script setup lang='ts'>
import {computed, ref, watch} from 'vue'
const split = '::';

const props = defineProps({
  label: String,
  hideEditBtn: {
    type: Boolean,
    default: false,
  },
  editable: {
    type: Boolean,
    default: false,
  },
  required: {
    type: Boolean,
    default: false,
  },
  modelValue: {
    type: String,
    default: "",
  },
  placeholder: {
    type: String,
    default: "",
  },
})

const localValue = ref(props.modelValue.split(split))

const inputType = computed(() => {
  return 'text'
})
const readonly = computed(() => {
  return !props.editable;
});

const active = ref(-1);

watch(() => active.value, (value) => {
  if (value === -1) {
    emit('update:modelValue', localValue.value.join(split));
  }
});

const emit = defineEmits(["update:modelValue"])

function addInputs() {
  localValue.value.push('');
}
function removeInputs(index: number) {
  localValue.value.splice(index, 1);
}
</script>

<style scoped lang='scss'>
.node-inputs-field {
  padding-left: 15px;
  .label-container {
    display: flex;
    align-items: center;
    .label {
      font-size: 14px;
      color: #6D6F84;
      margin-right: 5px;
    }
  }
  .value-container {
    display: flex;
    align-items: center;
    .value {
      &.active {
        border-color: #6D6F84;
      }
      background-color: #151728;
      border-radius: 6px;
      border-color: #151728;
      border-width: 2px;
      border-style: solid;
      box-sizing: border-box;
      height: 40px;
      line-height: 40px;
      color: #ececec;
      padding-right: 15px;
      margin-top: 6px;
      flex: 1;
      overflow: hidden;
      text-overflow: ellipsis;
      white-space: nowrap;
      display: flex;
      align-items: center;
      user-select: auto;
      input {
        width: 100%;
        flex: 1;
        outline: none;
        color: #ececec;
        margin-right: 10px;
        padding-left: 15px;
      }
    }
  }
}
</style>