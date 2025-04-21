<template>
    <v-snackbar location="top" text v-model="visible" :color="color">
        {{ msg }}
        <!-- 关闭按钮 -->
        <template v-slot:actions v-if="showClose">
            <v-btn color="#ffffff" variant="text" @click="close" size="large">
                <span class="mdi mdi-close-circle-outline"></span>
            </v-btn>
        </template>
    </v-snackbar>
</template>
<script setup lang="ts">
import { computed } from 'vue'
import { useSnackBarStore } from '@/stores/snackbar';
const snackBarStore = useSnackBarStore();

const visible = computed(() => {
    return snackBarStore.visible;
})
const showClose = computed(() => {
    return snackBarStore.showClose;
})
const msg = computed(() => {
    return snackBarStore.msg;
})
const color = computed(() => {
    return snackBarStore.color;
})
const close = () => {
    snackBarStore.closeSnackbar();
}
</script>