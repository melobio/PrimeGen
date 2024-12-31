import { ChatMessage } from '@azure/openai';

export interface AgentInterface {
  init(conversationUUID: string): void;
  send({
    userInput,
    description,
    optionInfo,
  }: {
    userInput: string;
    description?: string;
    optionInfo?: any;
  }): AsyncGenerator<ChatMessage, void, unknown>;
}
