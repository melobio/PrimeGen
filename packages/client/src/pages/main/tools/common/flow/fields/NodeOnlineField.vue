<template>
  <div class='node-online-field'>
    <div class='label-container'>
      <span class='label'>
        {{label}}
      </span>
    </div>
    <div class='value-container'>
      <div class='value' @mousedown.stop='() => {}'>
        <input
          :value="onlineState"
          :readonly='true'
          ref='inputRef'
          :class="{online: modelValue}"
        />
      </div>
      <v-btn
        :disabled="!enableRefresh"
        variant='text'
        :style="{color: '#6D6F84', opacity: 1}"
        icon='mdi-refresh'
        @click="refreshState"
      />
    </div>
  </div>
</template>

<script setup lang='ts'>
import {computed, ref} from 'vue'

const props = defineProps({
  label: String,
  modelValue: {
    type: Boolean,
    default: false,
  },
  enableRefresh: {
    type: Boolean,
    default: true,
  }
})

const inputRef = ref<HTMLInputElement>()
const onlineState = computed(() => {
  return props.modelValue ? 'online' : 'offline';
});

const emit = defineEmits(["update:modelValue", 'onRefresh'])
async function refreshState() {
  emit('onRefresh')
}

</script>

<style scoped lang='scss'>
.node-online-field {
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
      background-color: #151728;
      border-radius: 6px;
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
        color: #6D6F84;
        margin-right: 10px;
        padding-left: 15px;
        &.online {
          color: #ececec;
        }
      }
    }
  }
}
</style>