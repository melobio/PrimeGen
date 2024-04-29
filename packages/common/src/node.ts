// export const ToolTypeAgents = "Agents";
// export const ToolTypeRobot = "Devices";
// export const ToolTypeChains = "Chains";
// export const ToolTypeLoaders = "Loaders";
// export const ToolTypeEmbeddings = "Embeddings";
// export const ToolTypeInput = "Input";
// export const ToolTypeLLMs = "LLMs";
// export const ToolTypeMemories = "Memories";
// export const ToolTypeOutputParser = "OutputParser";
// export const ToolTypePrompts = "Prompts";
// export const ToolTypeRetrievers = "Retrievers";
// export const ToolTypeTextSplitters = "TextSplitters";
// export const ToolTypeToolsKits = "ToolsKits";
// export const ToolTypeTools = "Tools";
// export const ToolTypeUtilities = "Utilities";
// export const ToolTypeVectorStores = "VectorStores";
// export const ToolTypeWrappers = "Wrappers";
// import {NodeField} from "./node-field";

import type { NodeField } from "./node-field";
import type {FunctionItem} from "./functions";

export enum NodeTypes {
  Experiment = "Experiment",
  Agents = "Agents",
  Devices = "Devices",
  Chains = "Chains",
  Loaders = "Loaders",
  Embeddings = "Embeddings",
  Input = "Input",
  LLMs = "LLMs",
  Memories = "Memories",
  OutputParser = "OutputParser",
  Prompts = "Prompts",
  Retrievers = "Retrievers",
  TextSplitters = "TextSplitters",
  ToolsKits = "ToolsKits",
  Tools = "Tools",
  Utilities = "Utilities",
  VectorStores = "VectorStores",
  Wrappers = "Wrappers",
}

export const NodeTypeColors = {
  [NodeTypes.Experiment]: '#5E45B7',
  [NodeTypes.Agents]: '#8640B8',
  [NodeTypes.Devices]: '#4C66B9',
  [NodeTypes.Chains]: '#ED7D30',
  [NodeTypes.Loaders]: '#83AA51',
  [NodeTypes.Embeddings]: '#66B8A6',
  [NodeTypes.Input]: '#4BA3E3',
  [NodeTypes.LLMs]: '#5E45B7',
  [NodeTypes.Memories]: '#ECBA6A',
  [NodeTypes.OutputParser]: '#DCA945',
  [NodeTypes.Prompts]: '#4C66B9',
  [NodeTypes.Retrievers]: '#D6AE65',
  [NodeTypes.TextSplitters]: '#AC7EB2',
  [NodeTypes.ToolsKits]: '#C93E36',
  [NodeTypes.Tools]: '#3D485D',
  [NodeTypes.Utilities]: '#3D485D',
  [NodeTypes.VectorStores]: '#A5884C',
  [NodeTypes.Wrappers]: '#D33C79',
}

export interface Node {
  id: number;
  uuid: string;
  name: string;
  desc: string;
  shortDesc: string;
  type: NodeTypes;
  version: string;
  fields: NodeField[];
  functions: FunctionItem[];
  children?: NodeChild[];
  parent?: NodeChild;
}

export interface NodeChild {
  type: NodeTypes;
  name: string;
  uuid: string;
  required: boolean;
  functionName: string;
  functionDesc: string;
}