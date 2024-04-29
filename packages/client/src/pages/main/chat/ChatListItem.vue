<template>
  <div class='chat-list-item'>
    <div class='item-title' :class="{ active }">
      {{ name }}
    </div>
    <div class="del-icon" v-if="!disabledDelUUIDs.includes(uuid)" @click.stop="dialog = true">
      <v-icon icon="mdi-close-circle"></v-icon>
    </div>
    <div class='item-content'>
      {{ desc }}
    </div>
    <v-dialog v-model="dialog" width="auto">
      <v-card max-width="500"
        text="Are you sure you want to delete the experiment?" title="Delete Experiment">
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
import { computed, ref } from "vue";
const props = defineProps<Chat>()
const chatStore = useChatStore();
const active = computed(() => chatStore.currentChat?.uuid === props.uuid);
const dialog = ref(false);
const disabledDelUUIDs = computed(()=>{
  return [chatStore.chats[0].uuid,chatStore.chats[1].uuid]
})
const handleDeleteCon =async () => {
  await chatStore.handleDeleteConversationByUUID(props.uuid);
}
</script>

<style scoped lang='scss'>
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

}

.chat-list-item+.chat-list-item {
  margin-top: 10px;
}
</style>