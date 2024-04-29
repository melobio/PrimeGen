import {Node, NodeTypes} from "./node";

export enum PromptType {
  PromptTemplate = 'Prompt Template'
}
export enum PromptFields {
  PROMPT = 'prompt',
  ROLE = 'role',
}

export interface PromptNode extends Node {
  type: NodeTypes.Prompts;
  promptType: PromptType;
}