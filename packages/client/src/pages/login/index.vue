<template>
  <div class='login'>
    <LoginCard />
  </div>
</template>

<script setup lang='ts'>
import LoginCard from '@/pages/login/LoginCard.vue'
</script>

<style scoped lang='scss'>
.login {
  height: 100%;
  background: #151728;
  display: flex;
  align-items: center;
  justify-content: center;
}
</style>