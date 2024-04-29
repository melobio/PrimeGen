import { ChatMessage } from '@azure/openai';

export interface AgentInterface {
  init(conversationUUID: string): void;
  send(
    userInput: string,
    description: string,
  ): AsyncGenerator<ChatMessage, void, unknown>;
}
