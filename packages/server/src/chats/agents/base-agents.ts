import { AgentInterface } from './agent-interface';
import { Logger } from '@nestjs/common';
import {
  AzureKeyCredential,
  ChatCompletions,
  ChatMessage,
  OpenAIClient,
} from '@azure/openai';
import { OpenAI } from 'openai';
import { Agents } from '../../tools/agents/entities/agents.entity';
import { LLMFields, LLMType, NodeTypes } from '@xpcr/common';
import { LlmsEntity } from '../../tools/llms/entities/llms.entity';
import { DataSource } from 'typeorm';
import { Role } from '../entities/message.entity';
import { AgentMessageEntity } from '../entities/agent-message.entity';

export abstract class BaseAgents implements AgentInterface {
  get logger(): Logger {
    return this.getLogger();
  }
  needSummarize = true; // if need Router LLM to summarize
  client: OpenAIClient;
  GPTClient: OpenAI;
  llm_type: LLMType;
  conversationUUID: string;
  history: AgentMessageEntity[] = [];
  abstract getAvailableFunctions(): {
    [key: string]: (params: object) => AsyncGenerator<any, void, any>;
  };
  abstract getLogger(): Logger;
  abstract getInitMessages(): ChatMessage[];
  protected constructor(
    readonly agent: Agents,
    readonly dataSource: DataSource,
  ) {}
  protected *mockGenerating(text: string) {
    // mock generating
    for (let i = 0; i < text.length; i++) {
      yield {
        content: text[i],
      };
    }
  }

  async init(conversationUUID: string) {
    this.conversationUUID = conversationUUID;
    let apiKey = '';
    let apiBase = '';
    let apiEngine = '';
    const llmChildren = this.agent.children.filter(
      (nodeChild) => nodeChild.type === NodeTypes.LLMs,
    );
    // init azure llm
    for (const llmChild of llmChildren) {
      const llm = await this.dataSource.manager.findOne(LlmsEntity, {
        where: { uuid: llmChild.uuid },
      });
      if (llm && llm.llmType === LLMType.Azure) {
        llm.fields?.forEach((field) => {
          if (field.key === LLMFields.API_KEY) {
            apiKey = field.value;
          } else if (field.key === LLMFields.API_BASE) {
            apiBase = field.value;
          } else if (field.key === LLMFields.API_ENGINE) {
            apiEngine = field.value;
          }
        });
      }
    }
    if (!apiKey || !apiBase || !apiEngine) {
      apiKey = process.env.OPENAI_API_KEY;
      apiBase = process.env.OPENAI_API_BASE;
      apiEngine = process.env.OPENAI_API_ENGINE;
      this.logger.warn(`LLM config not found, use default config.`);
    }
    this.logger.log(`open ai key: ${apiKey}`);
    this.logger.log(`open ai base: ${apiBase}`);
    this.logger.log(`open ai engine: ${apiEngine}`);
    if (process.env.LLM_TYPE.toLowerCase() == 'azure') {
      this.llm_type = LLMType.Azure;
      this.client = new OpenAIClient(apiBase, new AzureKeyCredential(apiKey));
    } else if (process.env.LLM_TYPE.toLowerCase() == 'openai') {
      this.llm_type = LLMType.OpenAI;
      this.GPTClient = new OpenAI({
        baseURL: apiBase,
        apiKey: apiKey,
      });
    }
    this.logger.log(`init openai success.`);
    const agentMessages = await this.dataSource.manager.find(
      AgentMessageEntity,
      {
        where: {
          agentType: this.agent.agentType,
          conversationUUID,
        },
      },
    );
    this.history = agentMessages;
    this.logger.debug(
      `Agent#${this.agent.name} history(${agentMessages.length})`,
    );
  }

  async saveMessage(content: string, role: Role) {
    let agentInputMessage = new AgentMessageEntity();
    agentInputMessage.role = role;
    agentInputMessage.content = content;
    agentInputMessage.agentType = this.agent.agentType;
    agentInputMessage.conversationUUID = this.conversationUUID;
    agentInputMessage = await this.dataSource.manager.save(agentInputMessage);
    this.history.push(agentInputMessage);
    return agentInputMessage;
  }

  async *send({
    userInput,
    description,
    optionInfo,
  }: {
    userInput: string;
    description?: string;
    optionInfo?: any;
  }) {
    this.logger.log(`> userInput: ${JSON.stringify(userInput)}`);
    this.logger.log(`> description: ${description}`);
    const messages = [
      ...this.getInitMessages(),
      ...this.history.map<ChatMessage>((message) => {
        return {
          role: message.role === Role.None ? Role.Assistant : message.role,
          content: message.content,
        };
      }),
      { role: Role.User, content: userInput },
    ];
    // save user input
    await this.saveMessage(userInput, Role.User);
    const functions = this.agent.functions;
    let chatCompletions = null;
    if (this.llm_type == LLMType.Azure) {
      chatCompletions = await this.client.getChatCompletions(
        process.env.OPENAI_API_ENGINE,
        messages,
        {
          functionCall: 'auto',
          functions: functions,
        },
      );
    } else {
      const chatCompletionMessages: OpenAI.ChatCompletionMessageParam[] =
        messages.map((message) => {
          if ('name' in message) {
            return message;
          } else {
            return { ...message, name: '' }; // or some other default value
          }
        }) as OpenAI.ChatCompletionMessageParam[];
      const params: OpenAI.ChatCompletionCreateParams = {
        model: process.env.OPENAI_API_ENGINE,
        messages: chatCompletionMessages,
        stream: false,
      };
      chatCompletions = await this.GPTClient.chat.completions.create(params);
    }
    // this.logger.log(`> ${JSON.stringify(message)}`);
    // for await (const resp of this._handleChatCompletions(
    //   messages,
    //   chatCompletions,
    // )) {
    //   this.logger.log(`<< "${this.agent.name}" Response:`);
    //   console.warn(resp.content);
    //   yield resp;
    // }
    // yield* this._handleChatCompletions(messages, chatCompletions);
    let outputContent = '';
    for await (const msg of this._handleChatCompletions(
      messages,
      chatCompletions,
    )) {
      outputContent += msg?.content ?? '';
      yield msg;
    }
    await this.saveMessage(outputContent, Role.Assistant);
    // const resp = await this._handleChatCompletions(messages, chatCompletions);
    // this.logger.log(`< "${this.agent.name}" Response: ${resp.content}`);
    // return resp;
  }

  private async *_handleChatCompletions(
    messages: ChatMessage[],
    chatCompletions: ChatCompletions,
  ) {
    // parse OpenAI response ChatCompletions, call AvailableFunctions
    const responseMessage = chatCompletions.choices[0].message;
    if (responseMessage.functionCall) {
      const functionName = responseMessage.functionCall.name;
      this.logger.log(
        `> Use Function: ${functionName}, Arguments: ${JSON.stringify(
          responseMessage.functionCall.arguments,
        )}}`,
      );
      const availableFunctions = this.getAvailableFunctions();
      const functionToCall = availableFunctions[functionName];
      // this.logger.log(`functionToCall ${functionToCall}`);
      let functionResponse = '';
      // reply of agent
      for await (const msg of functionToCall(
        JSON.parse(responseMessage.functionCall.arguments),
      )) {
        if (msg.role === Role.Assistant) {
          // Agent message as AI response
          functionResponse = msg.content;
          continue;
        }
        yield msg;
      }
      // const functionResponse: { state: boolean } = await functionToCall(
      //   JSON.parse(responseMessage.functionCall.arguments),
      // );
      this.logger.log(
        `< ${functionName} return ${JSON.stringify(functionResponse)}`,
      );
      messages.push({
        role: responseMessage.role,
        functionCall: responseMessage.functionCall,
        content: undefined,
      });
      messages.push({
        role: 'function',
        name: functionName,
        content: functionResponse,
      });

      let chatCompletions = null;
      if (this.llm_type == LLMType.Azure) {
        chatCompletions = await this.client.getChatCompletions(
          process.env.OPENAI_API_ENGINE,
          messages,
          {
            functionCall: 'auto',
            functions: this.agent.functions,
          },
        );
      } else {
        const chatCompletionMessages: OpenAI.ChatCompletionMessageParam[] =
          messages.map((message) => {
            if ('name' in message) {
              return message;
            } else {
              return { ...message, name: '' }; // or some other default value
            }
          }) as OpenAI.ChatCompletionMessageParam[];
        const params: OpenAI.ChatCompletionCreateParams = {
          model: process.env.OPENAI_API_ENGINE,
          messages: chatCompletionMessages,
          stream: false,
        };
        chatCompletions = await this.GPTClient.chat.completions.create(params);
      }

      // yield await this._handleChatCompletions(messages, chatCompletions);
      for await (const msg of this._handleChatCompletions(
        messages,
        chatCompletions,
      )) {
        yield msg;
      }
    } else {
      yield chatCompletions.choices[0].message;
    }
  }

  // For query Agent
  protected async *sendQueryResult(searchResult: string, success: boolean) {
    if (success) {
      const content = {
        success: true,
        result: searchResult,
        retry: false,
      };
      yield {
        content: JSON.stringify(content),
        role: Role.Assistant,
      };
    } else {
      const content = {
        success: false,
        result: 'Query fail',
        retry: false,
      };
      yield {
        content: JSON.stringify(content),
        role: Role.Assistant,
      };
    }
  }
}
