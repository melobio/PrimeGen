<template>
  <div class='tool-list'>
    <div class='search-container'>
      <div class='search'>
        <v-icon icon='mdi-magnify' style='color: #545665'/>
        <input class='search-input' placeholder='Search for tools'/>
      </div>
    </div>

    <div class='list-container'>
      <ToolListItem
        v-for='(item) in toolsStore.allTools'
        :item='item'
        :active-item='toolsStore.currentTool'
        :key='item.name'
        @click='toolsStore.currentTool = item'
      />
    </div>
  </div>
</template>

<script setup lang='ts'>
import { useToolsStore } from '@/stores/tools'
import ToolListItem from '@/pages/main/tools/ToolListItem.vue'
const toolsStore = useToolsStore();
</script>

<style scoped lang='scss'>
.tool-list {
  width: 300px;
  height: 100%;
  padding: 15px 0;
  display: flex;
  flex-direction: column;
  .search-container {
    padding: 0 15px;
    margin-top: 25px;
    display: flex;
    flex-direction: row;
    .search {
      flex: 1;
      height: 40px;
      background: #1E2032;
      border-radius: 6px;
      display: flex;
      align-items: center;
      padding: 0 15px;
      .search-input {
        height: 100%;
        flex: 1;
        margin-left: 10px;
        border: none;
        outline: none;
        color: white;
        &::placeholder {
          color: #545665;
        }
      }
    }
  }

  .list-container {
    padding: 20px 15px 0;
    overflow-y: auto;
    flex: 1;
  }
}
</style>