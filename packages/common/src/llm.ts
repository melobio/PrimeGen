import {Node, NodeTypes} from "./node";

export enum LLMType {
  Azure = 'Azure',
  OpenAI = 'OpenAI',
}

export enum LLMFields {
  API_KEY = 'apiKey',
  API_BASE = 'apiBase',
  API_ENGINE = 'apiEngine',
}
export interface LLMNode extends Node {
  type: NodeTypes.LLMs;
  llmType: LLMType;
}