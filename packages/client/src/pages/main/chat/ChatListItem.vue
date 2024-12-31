<template>
  <div class='chat-list-item'>
    <div class='item-title' :class="{ active }">
      <span style="overflow: hidden;text-overflow: ellipsis;" v-if="!editTitle" :title="name">{{ name }}</span>
      <input v-show="editTitle" v-model="newTitle" ref="titleInput" type="text" @blur="handleEditTitle" @keydown.enter="handleEditTitle"
        class="title-input" />
    </div>
    <!-- <div class="del-icon" v-if="!disabledDelUUIDs.includes(uuid)" @click.stop="dialog = true">
      <v-icon icon="mdi-close-circle"></v-icon>
    </div> -->
    <div class="more-con" @click.stop>
      <v-menu transition="slide-y-transition">
        <template v-slot:activator="{ props }">
          <div class="more-icon" v-bind="props">
            <img width="15" height="15" src='@/assets/more-info.svg' alt='' />
          </div>
        </template>
        <v-list class="list-menu">
          <v-list-item value="b" @click="()=>{editTitle = true;newTitle=name}">
            <template v-slot:title>
              <div style="display:flex;align-items: center;">
                <v-icon style="font-size: 18px;margin-right: 5px;" icon="mdi mdi-pencil"></v-icon>
                <span>Rename</span>
              </div>
            </template>
          </v-list-item>
          <v-list-item value="a" @click="dialog = true">
            <template v-slot:title>
              <div style="display:flex;color:#FF2728;align-items: center;">
                <v-icon style="font-size: 18px;margin-right: 5px;"  icon="mdi-delete"></v-icon>
                <span>Delete</span>
              </div>
            </template>
          </v-list-item>
        </v-list>
      </v-menu>
    </div>
    <div class='item-content'>
      {{ desc }}
    </div>
    <v-dialog v-model="dialog" width="auto">
      <v-card max-width="500" text="Are you sure you want to delete the chat?" title="Delete Chat">
        <template v-slot:actions>
          <span class="ms-auto">
            <v-btn variant="text" color='#2db489' class="mr-2" @click='dialog = false'>
              Cancel
            </v-btn>
            <v-btn variant='elevated' color='red' @click='handleDeleteCon'>
              Confirm
            </v-btn>
          </span>
        </template>
      </v-card>
    </v-dialog>
  </div>
</template>

<script setup lang='ts'>
import type { Chat } from '@/stores/chat-types'
import { useChatStore } from "@/stores/chats";
import { computed, ref, watch, nextTick } from "vue";
import { useSnackBarStore } from '@/stores/snackbar';
const snackBarStore = useSnackBarStore();
const props = defineProps<Chat>()
const chatStore = useChatStore();
const active = computed(() => chatStore.currentChat?.uuid === props.uuid);
const dialog = ref(false);
const editTitle = ref(false);
const newTitle = ref('');
const titleInput = ref();

watch(() => editTitle.value, () => {
  if (editTitle.value) {
    console.log('titleInput.value===>', titleInput.value)
    nextTick(() => {
      titleInput.value.focus();
    })
  }
})
const disabledDelUUIDs = computed(() => {
  return [chatStore.chats[chatStore.chats.length - 2].uuid, chatStore.chats[chatStore.chats.length - 1].uuid]
})
const handleDeleteCon = async () => {
  await chatStore.handleDeleteConversationByUUID(props.uuid);
}
const handleEditTitle = async () => {
  if(newTitle.value){
    await chatStore.handleReNameConversationByUUID(props.uuid,newTitle.value);
  }else{
    snackBarStore.openSnackbar({
      msg: 'Chat name cannot be empty',
      color: 'red'
    })
  }
  editTitle.value = false;
}
</script>

<style scoped lang='scss'>
.list-menu{
  .list-item-icon{

  }
}
.chat-list-item {
  border-radius: 10px;
  background-color: #1E2032;
  //min-height: 100px;
  padding: 15px;
  transition: background-color 0.3s linear;
  cursor: pointer;
  display: flex;
  flex-direction: column;
  overflow: hidden;
  position: relative;

  .del-icon {
    opacity: 0;
    position: absolute;
    right: 4px;
    top: 5px;
    color: #545665;
    font-size: 12px;
    transition: all 0.3s ease-in-out;
  }

  &:hover {
    background-color: #292B3C;

    .del-icon {
      color: #616374;
      opacity: 1;
    }
  }

  .item-title {
    color: #ECECEC;
    font-size: 16px;
    text-overflow: ellipsis;
    overflow: hidden;
    white-space: nowrap;
    height: 23px;
    line-height: 23px;
    display: flex;
    align-items: center;

    &.active {
      color: #31DA9F;
    }
  }

  .item-content {
    margin-top: 5px;
    color: #545665;
    flex: 1;
    overflow: hidden;
    text-overflow: ellipsis;
    white-space: nowrap;
  }

  .more-con {
    position: absolute;
    height: 100%;
    width: 70px;
    right: 0;
    top: 0;
    opacity: 0;
    display: flex;
    justify-content: flex-end;
    align-items: center;
    padding: 15px;
    transition: all 0.3s ease-in-out;
    transform: translateX(10px);
    background: linear-gradient(to right, rgba(0, 0, 0, 0) 0%, #07070c 50%);

    .more-icon {
      padding: 5px;
      border-radius: 50%;
      display: flex;
      justify-content: center;
      align-items: center;
      transition: all 0.3s ease-in-out;

      &:hover {
        background-color: #545665;
      }
    }
  }

  &:hover {
    background-color: #292B3C;

    .more-con {
      transform: translateX(0px);
      opacity: .5;
    }
  }

  .title-input {
    width: 100%;
    border: 1px solid #2db489;
    padding: 0 5px;
    height: 100%;
    &:focus{outline:none;}
  }
}

.chat-list-item+.chat-list-item {
  margin-top: 10px;
}
</style>