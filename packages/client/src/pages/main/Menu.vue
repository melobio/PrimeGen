<template>
  <div class='menu'>
    <div class='user'>
      <img alt='' src='@/assets/user.png' class='icon'/>
      <span class='name'>
        X-Agents
      </span>
    </div>
    <div class='menu-list'>
      <template v-for='(sub, group) in menus' :key='group'>
        <template v-for='(menu, index) in sub' :key='menu.name'>
          <div class='menu-item-divider' :class='{ active: menu.name === modelValue }'/>
          <div class='menu-item-container' :class='{ active: menu.name === modelValue }'>
            <div class='menu-item'
                 @click='selectMenu(menu)'
                 :class='{ active: menu.name === modelValue }'>
              <tint-img class='menu-item-icon'
                        :src='menu.icon'
                        :tint="menu.name === modelValue ? '#31DA9F': 'white'"/>
              <span class='menu-item-name'>
                {{menu.name}}
              </span>
            </div>
          </div>
        </template>
        <div class='group-item-divider'/>
      </template>
    </div>
    <div class='logout-container'>
      <div class='flex-1-1'/>
      <v-btn
        color='#343648'
        elevation='0'
        :style="{color:'#FF2728'}"
        width='160px' height='40px'
        prepend-icon='mdi-exit-to-app'>
        <span style='text-transform: none'>
          Logout
        </span>
        <v-dialog
          v-model="logoutDialog"
          activator="parent"
          width="300"
        >
          <LogoutDialog @onCancel='logoutDialog = false'/>
        </v-dialog>
      </v-btn>
    </div>
  </div>
</template>

<script setup lang='ts'>
import * as _ from 'lodash'
import { ref } from 'vue'
import { Menus, MENUS } from '@/pages/main/menus'
import TintImg from '@/components/TintImg.vue'
import LogoutDialog from '@/pages/main/LogoutDialog.vue'

const props = defineProps({
  modelValue: {
    type: String,
    default: 'Chats'
  },
});
const emit = defineEmits(['update:modelValue'])

const menus = ref(
  _.groupBy(MENUS, 'group')
)
const logoutDialog = ref(false)

function selectMenu(menu: Menus) {
  emit('update:modelValue', menu.name)
}
</script>

<style scoped lang='scss'>
.menu {
  height: 100%;
  overflow-y: auto;
  display: flex;
  flex-direction: column;
  .user {
    align-self: stretch;
    display: flex;
    flex-direction: column;
    align-items: center;
    justify-content: center;
    padding-top: 40px;
    padding-bottom: 41px;
    background: #1E2032;
    .icon {
      width: 50px;
      height: 50px;
    }
    .name {
      margin-top: 15px;
      color: #ECECEC;
      height: 20px;
      line-height: 20px;
    }
  }
  .menu-list {
    padding-left: 0;
    .menu-item-container {
      padding-left: 20px;
      background: #1E2032;
      .menu-item {
        cursor: pointer;
        height: 40px;
        display: flex;
        flex-direction: row;
        align-items: center;
        padding-left: 25px;
        color: white;
        .menu-item-icon {
          margin-right: 15px;
        }
        .menu-item-name {
          font-size: 14px;
        }
        &.active {
          background: #151728;
          border-radius: 10px 0 0 10px;
          color: #31DA9F;
        }
      }
      &.active + div {
        border-top-right-radius: 10px;
      }
    }
    .menu-item-divider {
      height: 10px;
      background: #1E2032;
      &.active {
        border-bottom-right-radius: 10px;
      }
    }
    .group-item-divider {
      height: 30px;
      background: #1E2032;
    }
  }
  .logout-container {
    flex: 1;
    background: #1E2032;
    display: flex;
    flex-direction: column;
    align-items: center;
    padding: 20px 0;
  }
}
</style>