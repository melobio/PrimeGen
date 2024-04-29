<template>
  <v-dialog
    :model-value='modelValue'
    @close="$emit('update:modelValue', false)"
  >
    <div class='prompt-temp-edit-dialog'>
      <div class='title-container'>
        <span class='title'>{{title}}</span>
        <v-btn icon='mdi-close-circle' variant='plain'
               @click="$emit('update:modelValue', false)"
          color='#3F485B'>
        </v-btn>
      </div>
      <textarea
        class='edit-content'
        v-model='localTemp'
      />
      <div class='edit-footer'>
        <v-btn icon='mdi-alert-circle' variant='plain'
               color='#3F485B'>
        </v-btn>
        <div style='color: #3F485B; padding-top: 10px; flex: 1;'>
          <div>
            Prompt variable: {{}}
          </div>
          <div>
            The name of the prompt can be freely chosen within curly brackets, such as {variable_name}.
          </div>
        </div>
        <v-btn :color="'#3F485B'" :width='106' :height='40'
               class='align-self-center'
               @click="saveAndClose"
        >
          <span style='text-transform: none; color: #ECECEC; font-size: 16px;'>
            Save
          </span>
        </v-btn>
      </div>
    </div>
  </v-dialog>
</template>

<script setup lang='ts'>
import { onMounted, ref, watch } from 'vue'

const props = defineProps({
  title: String,
  modelValue: {
    type: Boolean,
    default: false,
  },
  template: {
    type: String,
    default: '',
  },
});

const localTemp = ref('')

watch([() => props.modelValue], ([modelValue]) => {
  if (modelValue) {
    localTemp.value = props.template;
  }
})

function onShow() {
  localTemp.value = props.template;
}

const emit = defineEmits([
  'update:modelValue',
  'update:template',
])
function saveAndClose() {
  emit('update:template', localTemp.value)
  emit('update:modelValue', false);
}
</script>

<style scoped lang='scss'>
.prompt-temp-edit-dialog {
  width: 1000px;
  height: 680px;
  background-color: #1E2032;
  border-radius: 10px;
  align-self: center;
  display: flex;
  flex-direction: column;
  .title-container {
    height: 63px;
    display: flex;
    justify-content: space-between;
    align-items: center;
    padding-left: 20px;
    padding-right: 10px;
    .title {
      font-size: 16px;
      color: #ECECEC;
    }
  }
  .edit-content {
    flex: 1;
    background-color: #151728;
    border-radius: 10px;
    color: #ECECEC;
    outline: none;
    padding: 20px;
    margin: 0 20px;
  }
  .edit-footer {
    height: 78px;
    display: flex;
    padding-right: 20px;
    padding-left: 5px;
  }
}
</style>