<template>
  <div class='node-select-field'>
    <div class='label-container'>
      <span class='label'>
        {{label}}
      </span>
      <v-icon icon='mdi-asterisk' size='12' color='red' v-if='required'/>
    </div>
    <div class='value-container'>
      <div class='value' @mousedown.stop='() => {}'>
        <span>{{localValue}}</span>
        <v-icon icon='mdi-unfold-more-horizontal'/>

        <v-menu activator="parent" offset='5'>
          <v-list style='background-color: #151728'>
            <v-list-item
              :active='item === localValue'
              :base-color="'#999999'"
              :color='"#ececec"'
              v-for="(item, index) in options"
              :key="item"
              :value="item"
              @click.once='saveItem(item)'
            >
              <v-list-item-title>{{ item }}</v-list-item-title>
            </v-list-item>
          </v-list>
        </v-menu>
      </div>
    </div>
  </div>
</template>

<script setup lang='ts'>
import { ref } from 'vue'
import type { PropType } from 'vue';

const props = defineProps({
  label: String,
  required: {
    type: Boolean,
    default: false,
  },
  modelValue: {
    type: String,
    default: "",
  },
  options: {
    type: Array as PropType<Array<string>>,
    default: () => ([]),
  }
})

const localValue = ref(props.modelValue)

const emit = defineEmits(["update:modelValue"])


function saveItem(item: string) {
  localValue.value = item;
  emit('update:modelValue', item);
}
</script>

<style scoped lang='scss'>
.node-select-field {
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
    cursor: pointer;
    .value {
      background-color: #151728;
      border-radius: 6px;
      height: 40px;
      line-height: 40px;
      color: #ececec;
      padding-left: 15px;
      padding-right: 10px;
      margin-top: 6px;
      flex: 1;
      overflow: hidden;
      text-overflow: ellipsis;
      white-space: nowrap;
      display: flex;
      justify-content: space-between;
      align-items: center;
      user-select: auto;
      input {
        width: 100%;
        outline: none;
        color: #ececec;
      }
    }
  }
}
</style>