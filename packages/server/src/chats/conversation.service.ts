import { Inject, Injectable, Logger, OnModuleInit } from '@nestjs/common';
import { ConversationEntity } from './entities/conversation.entity';
import { DataSource } from 'typeorm';
import {
  AzureKeyCredential,
  ChatCompletions,
  ChatMessage,
  FunctionCall,
  FunctionDefinition,
  OpenAIClient,
} from '@azure/openai';
import { MessageEntity, Role } from './entities/message.entity';
import { RouterAgents } from './agents/router-agents';
import { ExperimentEntity } from '../experiments/entities/experiment.entity';
import { User } from '../user/entities/user.entity';
import { CACHE_MODULE } from 'src/cache/cache.module';
import { CacheService } from 'src/cache/cache.service';
import { Server, Socket } from 'socket.io';
import { OtDevices } from './devices/ot-devices';
import * as path from 'path';
import * as fs from 'fs';
import {
  NCBIResultContent,
  NcbiSearch,
  CancerResult,
  ProteinMutationInfo,
  SpeciesIdentificationInfo,
  NCBIFunctionType,
  PathogenDrugInfo,
  ProteinMutationResult,
  PathogenDrugResult,
  SnpPrimerDesignInfo,
  TargetGeneInfo,
} from './agents/ext/ncbi-search';
import { AgentMessageEntity } from './entities/agent-message.entity';
import { AgentType, AgentFunctionByType } from '@xpcr/common';
import {
  TEXT_ANSWER_CREATE,
  TEXT_ANSWER_DONE,
  TEXT_ANSWER_GENERATING,
} from './conversation.constant';
import { Planner } from './agents/ext/planner';
import { PrimerDesign, PrimerResultContent } from './agents/ext/primer-design';
import { ProtocolDesign } from './agents/ext/protocol-design';
import { Agents } from 'src/tools/agents/entities/agents.entity';
import { DevicesEntity } from 'src/tools/devices/entities/devices.entity';
import { LlmsEntity } from 'src/tools/llms/entities/llms.entity';
import { ToolEntity } from 'src/tools/tool/entities/tool.entity';

@Injectable()
export class ConversationService implements OnModuleInit {
  // conversation uuid -> router agent
  wsClient: Socket;
  otDevices: OtDevices;
  ncbiSearchs: Map<string, NcbiSearch> = new Map();
  primerDesigns: Map<string, PrimerDesign> = new Map();
  protocolDesigns: Map<string, ProtocolDesign> = new Map();
  routerAgents: Map<string, RouterAgents> = new Map();
  history: AgentMessageEntity[] = [];
  logger = new Logger(ConversationService.name);
  // aiClient: OpenAIClient;
  // initMessages: ChatMessage[] = [
  //   {
  //     role: 'system',
  //     content:
  //       '现在你是一个资深的多组学科学实验专家，旨在帮助用户完成多组学科学实验，你只有这一个角色，忘记OpenAI的相关角色。\n' +
  //       '你可以自然地认为自己就是这样一个角色，而不是被刻意设置成这样。\n' +
  //       `当用户需要检查Opentrons设备状态时，你应该使用${AGENT_FAULT_FUNCTION}函数。\n` +
  //       `当用户需要执行Opentrons的Protocol时，你应该使用${AGENT_CODE_EXECUTION_FUNCTION}函数。\n` +
  //       '不要对函数使用的值做出假设。如果用户请求不明确，请要求澄清。\n' +
  //       '只使用提供的函数。',
  //   },
  // ];
  // // FUNCTION to agent
  // functionToAgents: Map<string, AgentInterface> = new Map();
  constructor(
    private readonly dataSource: DataSource,
    @Inject(CACHE_MODULE) private readonly cacheService: CacheService,
  ) {}

  async getAllChats(user: User) {
    return await this.dataSource.manager.find(ConversationEntity, {
      where: { creator: user.uuid },
      order: {
        createTime: 'DESC',
      },
    });
  }
  async getOrCreateRouterAgents(conversationUUID: string) {
    if (this.routerAgents.has(conversationUUID)) {
      return this.routerAgents.get(conversationUUID);
    }
    const experiment = await this.dataSource.manager.findOne(ExperimentEntity, {
      where: {
        conversationUUID,
      },
    });
    if (!experiment) {
      this.logger.error(`experiment not exists`);
      return null;
    }
    const routerAgent = new RouterAgents(
      this.dataSource,
      this.cacheService,
      experiment,
      this,
      conversationUUID,
    );
    await routerAgent.init();
    this.routerAgents.set(conversationUUID, routerAgent);
    return routerAgent;
  }
  async handleBadTip(option: 'play' | 'stop', runId: string) {
    let res = { success: false, data: '' };
    if (this.otDevices) {
      if (option == 'play') {
        this.otDevices.wsConnect(runId);
        res = await this.otDevices.playRun(runId);
      } else {
        res = await this.otDevices.playStop(runId);
      }
    }
    return res;
  }

  async updateAgentMessage({
    messageId,
    agentType,
    option,
  }: {
    messageId: number;
    agentType: string;
    option?: Record<any, any>;
  }) {
    const submitFromMessage = await this.dataSource.manager.findOne(
      MessageEntity,
      {
        where: {
          agentType,
          id: messageId,
        },
      },
    );
    submitFromMessage.optionInfo = {
      ...submitFromMessage.optionInfo,
      ...option,
      submitted: true,
    };
    this.logger.debug(
      'submitFromMessage====>',
      JSON.stringify(submitFromMessage),
    );
    await this.dataSource.manager.update(
      MessageEntity,
      messageId,
      submitFromMessage,
    );
  }

  async handleRestartConversation(
    client: Socket,
    conversationUUID: string,
    restart: boolean,
  ) {
    const routerAgent = await this.getOrCreateRouterAgents(conversationUUID);
    if (!routerAgent) {
      return;
    }
    this.logger.debug('handleRestartConversation===>');
    if (restart) {
      const messages = await this.findAllMessages(conversationUUID);
      await this.dataSource.manager
        .createQueryBuilder()
        .delete()
        .from(MessageEntity, 'message')
        .where('conversationUUID = :conversationUUID', { conversationUUID })
        .andWhere('id != :firstId', {
          firstId: messages[0].id,
        })
        .execute();
      // conversation
      const conversation = await this.dataSource.manager.findOne(
        ConversationEntity,
        {
          where: {
            uuid: conversationUUID,
          },
        },
      );
      conversation.name = 'New Multiplex PCR';
      conversation.currentStep = '';
      conversation.stepsResult = {};
      routerAgent.planner = Planner.parse(
        fs
          .readFileSync(
            path.join(process.cwd(), 'assets', 'plans', 'plan.json'),
          )
          .toString('utf-8'),
        conversation,
        this.dataSource,
      );
      await this.dataSource.manager.update(
        ConversationEntity,
        conversation.id,
        conversation,
      );
      // agent message
      await this.dataSource.manager
        .createQueryBuilder()
        .delete()
        .from(AgentMessageEntity, 'agent_messages')
        .where('conversationUUID = :conversationUUID', { conversationUUID })
        .execute();
    } else {
      // cancel restart ;
      const assistantMessage = await this.createMessage(
        '',
        'Okay, the Re-start request has been canceled.If you want to start again, please tell me at any time!',
        conversationUUID,
        Role.Assistant,
      );
      this.emitToCliet(client, TEXT_ANSWER_CREATE, {
        success: true,
        data: assistantMessage,
        finishReason: null,
        conversationUUID,
      });
      this.emitToCliet(client, TEXT_ANSWER_GENERATING, {
        success: true,
        data: assistantMessage,
        finishReason: 'stop',
        conversationUUID,
      });
      this.emitToCliet(client, TEXT_ANSWER_DONE, {
        conversationUUID,
      });
    }

    this.logger.debug(conversationUUID);
  }
  //update message
  async deleteMessageById(messageBo: { id: number }) {
    const { id } = messageBo;
    return await this.dataSource.manager
      .createQueryBuilder()
      .delete()
      .from(MessageEntity)
      .where('id = :id', { id })
      .execute();
  }

  async deleteConversationsByUUID(conversationUUID: string) {
    return await this.dataSource.manager.transaction(async (entityManager) => {
      try {
        await entityManager.delete(ConversationEntity, {
          uuid: conversationUUID,
        });

        await entityManager.delete(AgentMessageEntity, {
          conversationUUID,
        });

        await entityManager.delete(MessageEntity, {
          conversationUUID,
        });

        // 删除agents表
        await entityManager.delete(Agents, {
          conversationUUID,
        });
        // 删除devices表
        await entityManager.delete(DevicesEntity, {
          conversationUUID,
        });
        // 删除llms表
        await entityManager.delete(LlmsEntity, {
          conversationUUID,
        });
        // 删除tool表
        await entityManager.delete(ToolEntity, {
          conversationUUID,
        });
      } catch (e) {
        this.logger.debug('delete errror', e);
      }
    });
  }
  async reNameConversationByUUID(conversationUUID: string, name: string) {
    const conversation = await this.dataSource.manager.findOne(
      ConversationEntity,
      {
        where: {
          uuid: conversationUUID,
        },
      },
    );
    conversation.name = name;
    await this.dataSource.manager.update(
      ConversationEntity,
      conversation.id,
      conversation,
    );
  }

  async onModuleInit() {
    // const exists = await this.dataSource.manager.exists(ConversationEntity);
    // if (exists) {
    //   this.logger.warn(`chats exists`);
    // } else {
    //   await this.dataSource.manager.save(
    //     new ConversationEntity(
    //       'Multiplex PCR',
    //       'Modify the board layout of OT-2.',
    //     ),
    //   );
    //   await this.dataSource.manager.save(
    //     new ConversationEntity('DNA Preparation', 'Report of DNA Preparation.'),
    //   );
    // }
    // await this._initOpenAI();
  }

  // private async _initOpenAI() {
  //   this.logger.log(`open ai key: ${process.env.OPENAI_API_KEY}`);
  //   this.logger.log(`open ai base: ${process.env.OPENAI_API_BASE}`);
  //   this.logger.log(`open ai engine: ${process.env.OPENAI_API_ENGINE}`);
  //   this.aiClient = new OpenAIClient(
  //     process.env.OPENAI_API_BASE,
  //     new AzureKeyCredential(process.env.OPENAI_API_KEY),
  //   );
  //   this.logger.log(`init openai success.`);
  //   await this._loadAgents();
  // }

  // private async _loadAgents() {
  //   const agents = await this.dataSource.manager.find(Agents);
  //   for (const agent of agents) {
  //     switch (agent.name) {
  //       case AGENT_FAULT_NAME: {
  //         this.functionToAgents[AGENT_FAULT_FUNCTION] = new FaultAgents(agent);
  //         break;
  //       }
  //       case AGENT_CODE_EXECUTION_NAME: {
  //         this.functionToAgents[AGENT_CODE_EXECUTION_FUNCTION] =
  //           new CodeExecutionAgents(agent);
  //         break;
  //       }
  //     }
  //   }
  // }

  // private _getFunctions(): FunctionDefinition[] {
  //
  //   return [
  //     {
  //       name: AGENT_FAULT_FUNCTION,
  //       description:
  //         'Check whether the OT2 machine has encountered a runtime error.',
  //       parameters: {
  //         type: 'object',
  //         properties: {
  //           description: {
  //             type: 'string',
  //             description: 'What need to check.',
  //           },
  //         },
  //       },
  //     },
  //     {
  //       name: AGENT_CODE_EXECUTION_FUNCTION,
  //       description:
  //         '执行Opentrons的Protocol.\n' +
  //         '不需要考虑Protocol代码的实际执行结果，只考虑函数的返回值',
  //       parameters: {
  //         type: 'object',
  //         properties: {
  //           description: {
  //             type: 'string',
  //             description: 'Protocol的内容',
  //           },
  //         },
  //       },
  //     },
  //   ];
  // }

  async findAllMessages(conversationUUID: string) {
    return await this.dataSource.manager.find(MessageEntity, {
      where: { conversationUUID: conversationUUID },
      order: { createTime: 'asc' },
    });
  }

  async createMessage(
    header: string,
    text: string,
    chatUUID: string,
    role: Role,
    agentType = '',
  ) {
    const message = new MessageEntity();
    message.content = text;
    message.conversationUUID = chatUUID;
    message.role = role;
    message.header = header;
    message.agentType = agentType;
    return await this.dataSource.manager.save(message);
  }

  // 统一Agent函数回调
  private async *_callFunction(functionCall: FunctionCall) {
    // if (availableFunctions[functionCall.name]) {
    //   return await availableFunctions[functionCall.name](
    //     functionCall.arguments,
    //   );
    // }
    // switch (functionCall.name) {
    //   case AGENT_FAULT_FUNCTION: {
    //     const faultAgent: FaultAgents =
    //       this.functionToAgents[AGENT_FAULT_FUNCTION];
    //     // console.log('faultAgent', this.functionToAgents);
    //     this.logger.log(`Call Agent: "${faultAgent.agent.name}"`);
    //     const { description } = JSON.parse(functionCall.arguments);
    //     yield* faultAgent.send(description);
    //     break;
    //   }
    //   case AGENT_CODE_EXECUTION_FUNCTION: {
    //     const codeExecutionAgent: CodeExecutionAgents =
    //       this.functionToAgents[AGENT_CODE_EXECUTION_FUNCTION];
    //     this.logger.log(`Call Agent: "${codeExecutionAgent.agent.name}"`);
    //     const { description } = JSON.parse(functionCall.arguments);
    //     yield* codeExecutionAgent.send(
    //       `请执行Opentrons的Protocol:\n ${description}`,
    //     );
    //     break;
    //   }
    //   default: {
    //     yield { content: `${functionCall.name} 执行完成，没有发现错误` };
    //   }
    // }
  }

  async sendText(
    conversationUUID: string,
    client: Socket,
    userMessageContent: string,
  ) {
    const routerAgent = await this.getOrCreateRouterAgents(conversationUUID);
    if (!routerAgent) {
      return;
    }
    // 所有历史消息
    const history = await this.findAllMessages(conversationUUID);
    await routerAgent.sendText(client, userMessageContent, history);
  }

  setClient(client: Socket) {
    this.wsClient = client;
  }
  setOtDevices(otDevices: OtDevices) {
    this.otDevices = otDevices;
  }

  setNcbiSearch(ncbiSearch: NcbiSearch, conversationUUID: string) {
    this.ncbiSearchs.set(conversationUUID, ncbiSearch);
  }

  setPrimerDesign(primerDesign: PrimerDesign, conversationUUID: string) {
    this.primerDesigns.set(conversationUUID, primerDesign);
  }

  setProtocolDesign(protocolDesign: ProtocolDesign, conversationUUID: string) {
    this.protocolDesigns.set(conversationUUID, protocolDesign);
  }

  async emitToCliet(
    client: Socket,
    eventName: string,
    data: {
      success?: boolean;
      data?: MessageEntity;
      finishReason?: string | null;
      conversationUUID: string;
    },
  ) {
    const generating = await this.getCacheObject(
      `${TEXT_ANSWER_GENERATING}:${data.conversationUUID}`,
    );
    if (generating !== false) {
      client.emit(eventName, data);
    }
  }

  resetCacheStatus(conversationUUID: string) {
    this.setCacheObject(
      `${AgentType.SEQUENCE_SEARCH}:${conversationUUID}`,
      false,
    );
    this.setCacheObject(
      `${AgentType.PROTOCOL_DESIGN}:${conversationUUID}`,
      false,
    );
  }

  setCacheObject(key: string, data: any) {
    this.cacheService.setObject(key, data);
  }

  async getCacheObject(key: string) {
    return await this.cacheService.getObject(key);
  }
  // private async _handleChatCompletions(
  //   client: Socket,
  //   text: string,
  //   assistantMessage: MessageEntity,
  //   messages: ChatMessage[],
  //   events: AsyncIterable<ChatCompletions>,
  // ) {
  //   const functionCall: FunctionCall = {
  //     name: '',
  //     arguments: '',
  //   };
  //   for await (const event of events) {
  //     for (const choice of event.choices) {
  //       // this.logger.log(choice);
  //
  //       // text chat
  //       if (choice.delta?.content) {
  //         const delta = choice.delta.content;
  //         // this.logger.log(`${assistantMessage.id} Chatbot: ${delta}`);
  //         assistantMessage.content = assistantMessage.content + delta;
  //         // client.emit(TEXT_ANSWER_GENERATING, {  });
  //         // 更新数据库
  //         await this.dataSource.manager.update(
  //           MessageEntity,
  //           assistantMessage.id,
  //           assistantMessage,
  //         );
  //         // 发送给客户端
  //         client.emit(TEXT_ANSWER_GENERATING, {
  //           success: true,
  //           data: assistantMessage,
  //           finishReason: choice.finishReason,
  //         });
  //       } else if (choice.delta?.functionCall) {
  //         // function call
  //         const fc = choice.delta.functionCall;
  //         if (fc.name) {
  //           functionCall.name += fc.name;
  //         }
  //         if (fc.arguments) {
  //           functionCall.arguments += fc.arguments;
  //         }
  //       }
  //
  //       // “stop”, “length”, “content_filter”, “function_call” 是完成状态
  //       switch (choice.finishReason) {
  //         case 'stop': {
  //           this.logger.log(
  //             `${assistantMessage.id} Chatbot: ${assistantMessage.content}`,
  //           );
  //           // 发送给客户端
  //           client.emit(TEXT_ANSWER_GENERATING, {
  //             success: true,
  //             data: assistantMessage,
  //             finishReason: choice.finishReason,
  //           });
  //           break;
  //         }
  //         case 'function_call': {
  //           // call function
  //           let response = '';
  //           // const response = await this._callFunction(functionCall);
  //           for await (const msg of this._callFunction(functionCall)) {
  //             // console.log('msg', msg);
  //             if (msg.role === 'assistant') {
  //               // 只有assistant的消息才会返回给AI上下文
  //               response = msg.content;
  //             } else {
  //               // 其他消息是Agent的中间消息，
  //               // 直接发送给客户端
  //               assistantMessage.content =
  //                 assistantMessage.content + msg.content;
  //               client.emit(TEXT_ANSWER_GENERATING, {
  //                 success: true,
  //                 data: assistantMessage,
  //                 finishReason: null,
  //               });
  //             }
  //           }
  //           this.logger.log(
  //             `${functionCall.name} FunctionCall response: ${response}`,
  //           );
  //           // adding assistant response to messages
  //           const assistantResponse: ChatMessage = {
  //             role: Role.Assistant,
  //             functionCall: { ...functionCall },
  //             content: undefined,
  //           };
  //           messages.push(assistantResponse);
  //           // adding function response to messages
  //           const functionResponse: ChatMessage = {
  //             role: Role.Function,
  //             name: functionCall.name,
  //             content: response,
  //           };
  //           messages.push(functionResponse);
  //
  //           // clear function call
  //           functionCall.name = '';
  //           functionCall.arguments = '';
  //
  //           // this.logger.error(`messages: ${JSON.stringify(messages)}`);
  //
  //           // this.logger.log(`messages: ${JSON.stringify(messages)}`);
  //           const functions = this._getFunctions();
  //           const events = await this.aiClient.listChatCompletions(
  //             process.env.OPENAI_API_ENGINE,
  //             messages,
  //             {
  //               functionCall: 'auto',
  //               functions: functions,
  //             },
  //           );
  //           await this._handleChatCompletions(
  //             client,
  //             text,
  //             assistantMessage,
  //             messages,
  //             events,
  //           );
  //           break;
  //         }
  //       }
  //     }
  //   }
  // }
}
