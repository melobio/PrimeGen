<template>
  <div class='node-textarea-field'
    :style="{'padding-right': hideEditBtn ? '15px' : '0px'}"
  >
    <div class='label-container'>
      <span class='label'>
        {{label}}
      </span>
      <v-icon icon='mdi-asterisk' size='12' color='red' v-if='required'/>
    </div>
    <div class='value-container'>
      <div class='value' @mousedown.stop='() => {}'>
        <input
          v-model='localValue'
          :readonly='true'
          ref='inputRef'
          :placeholder="placeholder"
          :type='inputType'
        />
        <template v-if='secret'>
          <v-btn
            icon='mdi-eye-off'
            variant='text'
            :ripple='false'
            size='25'
            v-if='eyeOff'
            @click='eyeOff = false'
          />
          <v-btn
            icon='mdi-eye'
            variant='text'
            :ripple='false'
            size='25'
            v-else
            @click='eyeOff = true'
          />
        </template>
      </div>
      <template v-if="!hideEditBtn">
        <!--<v-btn-->
        <!--  v-if='editable && editing'-->
        <!--  variant='text'-->
        <!--  :style="{color: '#1484FC'}"-->
        <!--  icon='mdi-open-in-new'-->
        <!--  @click="exitEdit"-->
        <!--/>-->
        <v-btn
          variant='text'
          :style="{color: '#6D6F84'}"
          icon='mdi-open-in-new'
          @click="toEdit"
        />
      </template>
    </div>

    <PromptTempEditDialog
      v-model='templateEditDialog'
      v-model:template='localValue'
      :title="label"
    />
  </div>
</template>

<script setup lang='ts'>
import {computed, nextTick, ref, watch} from 'vue'
import PromptTempEditDialog from "@/pages/main/tools/common/flow/edit-dialog/PromptTempEditDialog.vue";

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
  secret: {
    type: Boolean,
    default: false,
  },
  placeholder: {
    type: String,
    default: "",
  },
})

const inputRef = ref<HTMLInputElement>()

const templateEditDialog = ref(false);

const localValue = ref(props.modelValue)
const eyeOff = ref(true) // can't see default

const inputType = computed(() => {
  return props.secret && eyeOff.value ? 'password' : 'text'
})

const emit = defineEmits(["update:modelValue"])

function toEdit() {
  // editing.value = true;
  templateEditDialog.value = true;
}

watch(() => templateEditDialog.value, (templateEditDialog) => {
  if (!templateEditDialog) {
    emit('update:modelValue', localValue.value);
  }
}, { deep: true })
</script>

<style scoped lang='scss'>
.node-textarea-field {
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
        color: #ececec;
        margin-right: 10px;
        padding-left: 15px;
      }
    }
  }
}
</style>