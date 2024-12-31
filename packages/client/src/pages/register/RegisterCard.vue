<template>
  <Buffer v-if="isDisplay" class="buffer-modal"/>
  <div class='register-card' :style="{'width': width + 'px', 'height': height + 'px'}">
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
        password placeholder='Enter password please'
        v-model='password'
      >
        <template v-slot:icon>
          <v-icon icon='mdi-lock'/>
        </template>
      </CustomInput>
    </div>
    <div class='input-container'>
      <CustomInput
        @onEnter='doRegister'
        password placeholder='Confirm password please'
        v-model='confirmPassword'
      >
        <template v-slot:icon>
          <v-icon icon='mdi-lock'/>
        </template>
      </CustomInput>
    </div>
    <div class="password-note">
      <ul v-if="errors.length > 0">
        <li v-for="error in errors" :key="error">{{ error }}</li>
      </ul>
    </div>
    <v-btn color='#292B3C' width='320px' height='44px' class='mt-10'
      @click='doRegister'
      :disabled="registerDisabled"
      >
      <span style='text-transform: none;'>
        Register
      </span>
    </v-btn>
    <v-btn color='#292B3C' width='320px' height='44px' class='mt-10'
      @click='toLoginPage'>
      <span style='text-transform: none;'>
        Cancel
      </span>
    </v-btn>
    <div class='flex-1-1'/>
    <div class='copy-right'>
      Copyright © 2024 MGI Inc.
    </div>
  </div>
</template>

<script setup lang='ts'>
import CustomInput from '@/pages/login/CustomInput.vue'
import Buffer from '@/pages/register/Buffer.vue'
import { computed, ref } from 'vue'
import { login,register } from '@/api'
import { useRouter } from 'vue-router'
const router = useRouter()
const baseURL = import.meta.env.BASE_URL;

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
const confirmPassword = ref('')
const isDisplay = ref(false)

const errors = computed(() => {
  const errors = [];

  if (!userName.value) {
    errors.push('Username cannot be empty');
  }

  // Check if password is at least 6 characters long
  if (password.value.length < 6) {
    errors.push('Password must be at least 6 characters long');
  }

  // Check if password contains at least one uppercase letter
  if (!/[A-Z]/.test(password.value)) {
    errors.push('Password must contain at least one uppercase letter');
  }

  // Check if password contains at least one lowercase letter
  if (!/[a-z]/.test(password.value)) {
    errors.push('Password must contain at least one lowercase letter');
  }

  // Check if password contains at least one number
  if (!/\d/.test(password.value)) {
    errors.push('Password must contain at least one number');
  }

  // Check if password and confirm password match
  if (password.value !== confirmPassword.value) {
    errors.push('Passwords do not match');
  }

  return errors;
})

const registerDisabled = computed(() => {
  return errors.value.length > 0
})

async function doRegister() {
  const { success } = await register<{ token: string }>({
    loginType: "userName",
    userName: userName.value,
    password: password.value,
    confirmPassword:confirmPassword.value
  });
  if (success) {
    isDisplay.value = true
    setTimeout(() => {
      // window redirect to Login page
      window.location.href = `${baseURL}/login`;
    }, 2000);
  }
}
async function toLoginPage(){
    await router.push({name:'Login'})
}
</script>

<style scoped lang='scss'>
.buffer-modal {
  position: fixed;
  z-index: 1000; 
}
.register-card {
  background: #1E2032;
  border-radius: 10px;
  display: flex;
  flex-direction: column;
  align-items: center;
  padding-top: 60px;
  .user-icon {
    width: 50px;
    height: 60px;
  }
  .primegen {
    color: #ECECEC;
    margin-top: 15px;
    margin-bottom: 60px;
  }
  .input-container {
    width: 320px;
  }
  .password-note {
    white-space: pre;
    display: block;
    align-content: start;
    width: 320px;
    font-style: italic;
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