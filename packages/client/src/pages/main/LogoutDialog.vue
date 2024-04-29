<template>
  <div class='login-dialog'>
    <div class='title'>
      Tips
    </div>
    <div class='content'>
      Are you sure to log out？
    </div>
    <div class='text-right mt-4'>
      <v-btn variant='text' color='#6D6F84' @click="$emit('onCancel')">
        <span style='text-transform: none'>
          Cancel
        </span>
      </v-btn>
      <v-btn variant='text' class='ml-3' color='red' @click='exit'>
        <span style='text-transform: none'>
          Confirm Exit
        </span>
      </v-btn>
    </div>
  </div>
</template>

<script setup lang='ts'>
import Cookies from 'js-cookie'
import { useRouter } from 'vue-router'
import {useChatStore} from "@/stores/chats";
const router = useRouter()
const { logout } = useChatStore();

defineEmits(['onCancel'])
function exit() {
  Cookies.remove('x-token');
  logout();
  router.push({ name: 'Login' });
}
</script>

<style scoped lang='scss'>
.login-dialog {
  background-color: #1E2032;
  border-radius: 10px;
  padding: 10px 20px 10px;
  .title {
    color: #ECECEC;
    font-size: 18px;
    line-height: 44px;
    height: 44px;
  }
  .content {
    color: #ECECEC;
    font-size: 15px;
    line-height: 44px;
    height: 44px;
  }
}
</style>