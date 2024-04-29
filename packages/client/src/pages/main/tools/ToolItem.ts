import Agents from '@/assets/tools/agents.png';
import Robot from '@/assets/tools/robot.png';
import Chains from '@/assets/tools/chains.png';
import Loaders from '@/assets/tools/loaders.png';
import Embeddings from '@/assets/tools/embeddings.png';
import Input from '@/assets/tools/input.png';
import LLMs from '@/assets/tools/llms.png';
import Memories from '@/assets/tools/memories.png';
import OutputParser from '@/assets/tools/outputparser.png';
import Prompts from '@/assets/tools/prompts.png';
import Retrievers from '@/assets/tools/retrievers.png';
import TextSplitters from '@/assets/tools/textsplitters.png';
import ToolsKits from '@/assets/tools/toolskits.png';
import Tool from '@/assets/tools/tools.png';
import Utilities from '@/assets/tools/utilities.png';
import VectorStores from '@/assets/tools/vectorstores.png';
import Wrappers from '@/assets/tools/wrappers.png';
import { NodeTypes } from '@xpcr/common/src/node';
import {NodeTypeColors} from "@xpcr/common/src/node";

export type ToolsName = NodeTypes

export class ToolItem {
  color: string;
  constructor(
    readonly name: NodeTypes,
    readonly icon: string,
  ) {
    this.color = NodeTypeColors[name];
  }
}

export const TOOLS: ToolItem[] = [
  new ToolItem(NodeTypes.Agents, Agents),
  new ToolItem(NodeTypes.Devices, Robot),
  new ToolItem(NodeTypes.Chains, Chains),
  new ToolItem(NodeTypes.Loaders, Loaders),
  new ToolItem(NodeTypes.Embeddings, Embeddings),
  new ToolItem(NodeTypes.Input, Input),
  new ToolItem(NodeTypes.LLMs, LLMs),
  new ToolItem(NodeTypes.Memories, Memories),
  new ToolItem(NodeTypes.OutputParser, OutputParser),
  new ToolItem(NodeTypes.Prompts, Prompts),
  new ToolItem(NodeTypes.Retrievers, Retrievers),
  new ToolItem(NodeTypes.TextSplitters, TextSplitters),
  new ToolItem(NodeTypes.ToolsKits, ToolsKits),
  new ToolItem(NodeTypes.Tools, Tool),
  new ToolItem(NodeTypes.Utilities, Utilities),
  new ToolItem(NodeTypes.VectorStores, VectorStores),
  new ToolItem(NodeTypes.Wrappers, Wrappers),
]