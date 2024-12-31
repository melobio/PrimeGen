import { ref, computed } from 'vue'
import { defineStore } from 'pinia'

export const useSnackBarStore = defineStore('snackbar', () => {
    const msg = ref('')
    const color = ref('')
    const visible = ref(false)
    const showClose = ref(true)
    const timeout = ref(5000)

    function closeSnackbar() {
        visible.value = false
    }
    function setShowClose(isShow: boolean) {
        showClose.value = isShow
    }
    function openSnackbar(options: { msg: string, color: string }) {
        visible.value = true
        msg.value = options.msg
        color.value = options.color
        setTimeout(() => {
            closeSnackbar();
        }, timeout.value);
    }
    return {
        msg,
        color,
        visible,
        showClose,
        timeout,
        closeSnackbar,
        setShowClose,
        openSnackbar
    }
})
