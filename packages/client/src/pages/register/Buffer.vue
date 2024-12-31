<template>
    <div class='login-dialog'>
        <img class="success-svg" :src='checkMark' alt='check mark'/>
        <br>
        <div class='content'>
            <p>Successfully registered.</p>
            <p>Redirecting to login......</p>
        </div>

    </div>
</template>

<script lang="ts" setup>
import checkMark from '@/assets/success-green-check-mark-icon.svg';
</script>

<style scoped lang='scss'>
.login-dialog {
    width: 400px;
    height: 600px;
    background-color: #1E2032;
    border-radius: 10px;
    padding: 10px 20px 10px;
    display: grid;
    align-items: center;
    align-content: center;
    justify-content: space-evenly;
    justify-items: center;

    .success-svg {
        width: 100px;
        height: 100px;
        display: block;
    }

    .content {
        color: #ECECEC;
        font-size: 15px;
        line-height: 44px;
        height: 44px;
        white-space: pre-line;
    }
}
</style>