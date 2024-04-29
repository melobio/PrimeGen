import { markRaw } from 'vue'
import PromptNodeComponent from '@/pages/main/tools/common/flow/PromptNode.vue'
import type { ToolsName } from '@/pages/main/tools/ToolItem'
import CommonNode from '@/pages/main/tools/common/flow/CommonNode.vue'
import type { NodeComponent, NodeTypesObject } from '@vue-flow/core'
import AgentNodeComponent from "@/pages/main/tools/common/flow/AgentNode.vue";
import LLMsNodeComponent from "@/pages/main/tools/common/flow/LLMNode.vue";
import DeviceNodeComponent from "@/pages/main/tools/common/flow/DeviceNode.vue";
import InputNodeComponent from "@/pages/main/tools/common/flow/InputNode.vue";
import ExperimentNodeComponent from "@/pages/main/tools/common/flow/ExperimentNode.vue";
import ToolNodeComponent from "@/pages/main/tools/common/flow/ToolNode.vue";
import { NodeTypes as NT} from '@xpcr/common/src/node';
import type {
  AgentNode,
  LLMNode,
  DeviceNode,
  PromptNode,
  InputNode,
  Node,
  ExperimentInterface,
  ToolNode
} from '@xpcr/common';
export const NodeTypes: NodeTypesObject = {
  [NT.Experiment]: markRaw(ExperimentNodeComponent) as NodeComponent,
  [NT.Agents]: markRaw(AgentNodeComponent) as NodeComponent,
  [NT.Devices]: markRaw(DeviceNodeComponent) as NodeComponent,
  [NT.Chains]: markRaw(CommonNode) as NodeComponent,
  [NT.Loaders]: markRaw(CommonNode) as NodeComponent,
  [NT.Embeddings]: markRaw(CommonNode) as NodeComponent,
  [NT.Input]: markRaw(InputNodeComponent) as NodeComponent,
  [NT.LLMs]: markRaw(LLMsNodeComponent) as NodeComponent,
  [NT.Memories]: markRaw(CommonNode) as NodeComponent,
  [NT.OutputParser]: markRaw(CommonNode) as NodeComponent,
  [NT.Prompts]: markRaw(PromptNodeComponent) as NodeComponent,
  [NT.Retrievers]: markRaw(CommonNode) as NodeComponent,
  [NT.TextSplitters]: markRaw(CommonNode) as NodeComponent,
  [NT.ToolsKits]: markRaw(CommonNode) as NodeComponent,
  [NT.Tools]: markRaw(ToolNodeComponent) as NodeComponent,
  [NT.Utilities]: markRaw(CommonNode) as NodeComponent,
  [NT.VectorStores]: markRaw(CommonNode) as NodeComponent,
  [NT.Wrappers]: markRaw(CommonNode) as NodeComponent,
}

export interface BaseProp {
  id: string;
  type: ToolsName;
  position: { x: number; y: number; };
  data: any;
}

export interface ExperimentProp extends BaseProp{
  type: ToolsName;
  data: ExperimentInterface;
}

export interface CommonProp extends BaseProp{
  type: ToolsName;
  data: Node;
}

export interface PromptProp extends BaseProp{
  type: NT.Prompts;
  data: PromptNode;
}

export interface AgentProp extends BaseProp{
  type: NT.Agents;
  data: AgentNode;
}
export interface LLMProp extends BaseProp{
  type: NT.LLMs;
  data: LLMNode;
}
export interface DeviceProp extends BaseProp{
  type: NT.Devices;
  data: DeviceNode;
}

export interface InputProp extends BaseProp{
  type: NT.Input;
  data: InputNode;
}

export interface ToolProp extends BaseProp{
  type: NT.Tools;
  data: ToolNode;
}

export type AllNodeProp = BaseProp | PromptProp | AgentProp | DeviceProp | InputProp;
