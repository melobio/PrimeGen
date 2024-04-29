import type {Ref} from 'vue'
import {computed, ref} from 'vue'
import {defineStore} from 'pinia'
import type {ToolItem, ToolsName} from '@/pages/main/tools/ToolItem'
import {TOOLS} from '@/pages/main/tools/ToolItem'
import {
  toolsAgentTemplates,
  toolsDeviceTemplates,
  toolsInputTemplates,
  toolsLLMTemplates,
  toolsPromptTemplates, toolsToolTemplates
} from '@/api'
import type {Node} from '@xpcr/common/src/node';
import {NodeTypes} from '@xpcr/common/src/node';
import type {AgentNode} from '@xpcr/common/src/agent';
import type {DeviceNode, InputNode, LLMNode, PromptNode, ToolNode} from "@xpcr/common";

type ToolsState = Map<ToolsName, Ref<Node[]>>;

export const useToolsStore = defineStore('tools', () => {
  const allTools = ref(TOOLS);
  const currentTool = ref<ToolItem>(TOOLS[0]);
  // console.log('NodeTypes', Object.values(NodeTypes))

  const toolsState: ToolsState = new Map([
    [NodeTypes.Agents, ref<AgentNode[]>([])],
    [NodeTypes.Devices, ref<Node[]>([])],
    [NodeTypes.Chains, ref<Node[]>([])],
    [NodeTypes.Loaders, ref<Node[]>([])],
    [NodeTypes.Embeddings, ref<Node[]>([])],
    [NodeTypes.Input, ref<Node[]>([])],
    [NodeTypes.LLMs, ref<Node[]>([])],
    [NodeTypes.Memories, ref<Node[]>([])],
    [NodeTypes.OutputParser, ref<Node[]>([])],
    [NodeTypes.Prompts, ref<Node[]>([])],
    [NodeTypes.Retrievers, ref<Node[]>([])],
    [NodeTypes.TextSplitters, ref<Node[]>([])],
    [NodeTypes.ToolsKits, ref<Node[]>([])],
    [NodeTypes.Tools, ref<Node[]>([])],
    [NodeTypes.Utilities, ref<Node[]>([])],
    [NodeTypes.VectorStores, ref<Node[]>([])],
    [NodeTypes.Wrappers, ref<Node[]>([])],
  ]);

  async function getCurrentToolList() {
    const toolList = toolsState.get(currentTool.value.name)!!;
    switch (currentTool.value.name) {
      case NodeTypes.Agents: {
        const { success, data } =  await toolsAgentTemplates<AgentNode[]>();
        toolList.value = data || [];
        break;
      }
      case NodeTypes.Devices: {
        const { success, data } =  await toolsDeviceTemplates<DeviceNode[]>();
        toolList.value = data || [];
        break;
      }
      case NodeTypes.Input: {
        const { success, data } =  await toolsInputTemplates<InputNode[]>();
        toolList.value = data || [];
        break;
      }
      case NodeTypes.LLMs: {
        const { success, data } =  await toolsLLMTemplates<LLMNode[]>();
        toolList.value = data || [];
        break;
      }
      case NodeTypes.Prompts: {
        const { success, data } =  await toolsPromptTemplates<PromptNode[]>()
        toolList.value = data || [];
        break;
      }
      case NodeTypes.Tools: {
        const { success, data } =  await toolsToolTemplates<ToolNode[]>()
        toolList.value = data || [];
        break;
      }
      default:
        break;
    }
  }

  const currentToolList = computed(() => {
    return toolsState.get(currentTool.value.name) || ref<Node[]>([]);
  })

  return {
    allTools,
    currentTool,
    getCurrentToolList,
    // subset
    toolsState,
    currentToolList,
  };
})