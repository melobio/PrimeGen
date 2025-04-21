<template>
  <div class='login-card' :style="{'width': width + 'px', 'height': height + 'px'}">
    <img class='user-icon' src='@/assets/user.png' alt=''>
    <span class='x-agents'>
      X-Agents
    </span>
    <div class='input-container'>
      <CustomInput
        placeholder='Enter account please'
        v-model='userName'
      >
        <template v-slot:icon>
          <v-icon icon='mdi-account-outline'/>
        </template>
      </CustomInput>
    </div>
    <div class='input-container'>
      <CustomInput
        @onEnter='doLogin'
        password placeholder='Enter password please'
        v-model='password'
      >
        <template v-slot:icon>
          <v-icon icon='mdi-lock'/>
        </template>
      </CustomInput>
    </div>
    <v-btn color='#292B3C' width='320px' height='44px' class='mt-11'
      @click='doLogin'>
      <span style='text-transform: none;'>
        Login
      </span>
    </v-btn>
    <div class='flex-1-1'/>
    <div class='copy-right'>
      Copyright © 2023 MGI Inc.
    </div>
  </div>
</template>

<script setup lang='ts'>
import CustomInput from '@/pages/login/CustomInput.vue'
import { ref } from 'vue'
import { login } from '@/api'
import { useRouter } from 'vue-router'
const router = useRouter()
const props = defineProps({
  width: {
    type: Number,
    default: 400,
  },
  height: {
    type: Number,
    default: 600,
  }
})
const userName = ref('')
const password = ref('')
async function doLogin() {
  const { success } = await login<{ token: string }>({
    loginType: "userName",
    userName: userName.value,
    password: password.value,
  });
  if (success) {
    console.log('Login success.')
    await router.push({ name: 'Home' })
  }
}
</script>

<style scoped lang='scss'>
.login-card {
  background: #1E2032;
  border-radius: 10px;
  display: flex;
  flex-direction: column;
  align-items: center;
  padding-top: 60px;
  .user-icon {
    width: 50px;
    height: 50px;
  }
  .x-agents {
    color: #ECECEC;
    margin-top: 15px;
    margin-bottom: 60px;
  }
  .input-container {
    width: 320px;
  }
  .input-container + .input-container {
    margin-top: 30px;
  }
  .copy-right {
    color: #4F5159;
    margin-bottom: 15px;
  }
}
</style>