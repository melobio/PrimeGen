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
} from './agents/ext/ncbi-search';
import { AgentMessageEntity } from './entities/agent-message.entity';
import { AgentType, AgentFunctionByType } from '@xpcr/common';
import {
  TEXT_ANSWER_CREATE,
  TEXT_ANSWER_DONE,
  TEXT_ANSWER_GENERATING,
} from './conversation.constant';

@Injectable()
export class ConversationService implements OnModuleInit {
  // conversation uuid -> router agent
  wsClient: Socket;
  otDevices: OtDevices;
  ncbiSearch: NcbiSearch;
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
    return await this.dataSource.manager.findBy(ConversationEntity, {
      creator: user.uuid,
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

  objectToMarkdownTable = (data) => {
    const keys = Object.keys(data);
    const rows = Object.keys(data[keys[0]]).map((key) => {
      return keys.map((k) => data[k][key]);
    });

    const header = `| ${keys.join(' | ')} |\n`;
    const divider = `|${keys.map(() => ' :---: ').join('|')}|\n`;
    const body = rows.map((row) => `| ${row.join(' | ')} |`).join('\n');

    return `${header}${divider}${body}`;
  };

  checkIfAllEmpty = (obj) => {
    for (const key in obj) {
      if (Object.keys(obj[key]).length !== 0) {
        return false;
      }
    }
    return true;
  };

  async handleNcbiSearch({
    option,
    conversationUUID,
  }: {
    option: SpeciesIdentificationInfo | ProteinMutationInfo | PathogenDrugInfo;
    conversationUUID: string;
  }) {
    let res = { success: false, data: null };
    const routerAgent = await this.getOrCreateRouterAgents(conversationUUID);
    if (!routerAgent) {
      return;
    }
    if (this.ncbiSearch) {
      // 创建一个search agent的消息返回给client
      const agentMessage = await this.createMessage(
        '',
        '',
        conversationUUID,
        Role.Assistant,
        AgentType.SEQUENCE_SEARCH,
      );
      this.wsClient.emit(TEXT_ANSWER_CREATE, {
        success: true,
        data: agentMessage,
        finishReason: null,
      });
      this.wsClient.emit(TEXT_ANSWER_GENERATING, {
        success: true,
        data: agentMessage,
        finishReason: null,
        conversationUUID,
      });
      const agentMessages = await this.dataSource.manager.find(
        AgentMessageEntity,
        {
          where: {
            agentType: AgentType.SEQUENCE_SEARCH,
            conversationUUID,
          },
        },
      );
      const chatHistory = await this.findAllMessages(conversationUUID);
      const chatHistoryArr = chatHistory.map<ChatMessage>((item) => {
        return { role: item.role, content: item.content };
      });
      const history: AgentMessageEntity[] = agentMessages;
      const operationsObj = {};
      if (option?.operations) {
        option.operations.forEach((item) => {
          operationsObj[item.key] = item.value;
        });
      }
      const sumitOption = option
        ? { ...option, ...operationsObj, conversation: chatHistoryArr }
        : { ...operationsObj, conversation: chatHistoryArr };
      const nCBIResult: NCBIResultContent =
        await this.ncbiSearch.callFromService(
          JSON.stringify(sumitOption),
          history,
        );
      if (typeof nCBIResult == 'string') {
        agentMessage.content = `${nCBIResult} \n`;
        res = { success: true, data: nCBIResult };
      } else {
        const responses = nCBIResult.responses;
        agentMessage.content = `${responses.response} \n`;
        try {
          if (
            responses?.search_type ==
            NCBIFunctionType.species_identification_type
          ) {
            const non_target_path = responses?.non_target_path || [];
            const fileList = non_target_path.map((p) => {
              const url = path.join('/v1/files/xsearchdev', p);
              return `[${path.basename(url)}](${url})`;
            });
            agentMessage.content += `\n ${fileList.join('\n')}\n`;
          } else if (
            responses?.search_type ==
            NCBIFunctionType.pathogen_drug_resistance_type
          ) {
            const responses: PathogenDrugResult = nCBIResult.responses;
            const filter_result_path = responses?.filter_result_path || '';
            if (filter_result_path) {
              const path_url = path.join(
                '/v1/files/xsearchdev',
                filter_result_path,
              );
              const path_url_mdStr = `[${path.basename(
                path_url,
              )}](${path_url})`;
              agentMessage.content += `\n ${path_url_mdStr} ;`;
            }
            const sequences = responses?.sequences || '';
            if (sequences) {
              const path_url = path.join('/v1/files/xsearchdev', sequences);
              const path_url_mdStr = `[${path.basename(
                path_url,
              )}](${path_url})`;
              agentMessage.content += `\n ${path_url_mdStr} ;`;
            }
          } else if (
            responses?.search_type == NCBIFunctionType.protein_mutation_type
          ) {
            const responses: ProteinMutationResult = nCBIResult.responses;
            const protein_gene_seq_path =
              responses?.protein_gene_seq_path || [];
            const fileList = protein_gene_seq_path.map((p) => {
              const url = path.join('/v1/files/xsearchdev', p);
              return `[${path.basename(url)}](${url})`;
            });
            agentMessage.content += `\n ${fileList.join('\n')}\n`;
          } else if (responses?.search_type == NCBIFunctionType.cancer_type) {
            const responses: CancerResult = nCBIResult.responses;
            const target_gene_info = responses?.target_gene_info;
            if (Array.isArray(target_gene_info)) {
              const fileList = target_gene_info.map((p) => {
                const url = path.join('/v1/files/xsearchdev', p);
                return `[${path.basename(url)}](${url})`;
              });
              agentMessage.content += `\n ${fileList.join('\n')}\n`;
            } else if (target_gene_info) {
              const tableStr = this.objectToMarkdownTable(target_gene_info);
              agentMessage.content += tableStr;
            }
          }
        } catch (error) {
          this.logger.debug('handleNcbiSearch======>error');
          this.logger.debug(error);
        }
        agentMessage.optionInfo = responses;
        res = { success: true, data: responses };
      }
      // if (responses?.search_type == NCBIFunctionType.genetic_disorder_type) {
      //   const responses: TargetGeneInfo = nCBIResult.responses;
      //   const target_gene_info = responses.target_gene_info;
      //   if (this.checkIfAllEmpty(target_gene_info)) {
      //     agentMessage.content +=
      //       'The analysis results for the related genes were not found.';
      //   } else {
      //     const preDes = `${
      //       responses.response ||
      //       'these are the analysis results of the related genes.'
      //     }:\n `;
      //     const tableStr =
      //       preDes + this.objectToMarkdownTable(target_gene_info);
      //     agentMessage.content += tableStr;
      //   }
      // } else if (
      //   responses?.search_type == NCBIFunctionType.download_type
      // ) {
      //   const responses: TargetProteinMutationInfo = nCBIResult.responses;
      //   const protein_gene_seq_path = responses?.protein_gene_seq_path || [];
      //   const fileList = protein_gene_seq_path.map((p) => {
      //     const url = path.join('/v1/files/xsearchdev', p);
      //     return `[${path.basename(url)}](${url})`;
      //   });
      //   agentMessage.content = `${responses.response}\n${fileList.join(
      //     '\n',
      //   )}\n`;
      // }

      this.wsClient.emit(TEXT_ANSWER_GENERATING, {
        success: true,
        data: agentMessage,
        finishReason: null,
        conversationUUID,
      });

      await this.dataSource.manager.update(
        MessageEntity,
        agentMessage.id,
        agentMessage,
      );

      this.wsClient.emit(TEXT_ANSWER_GENERATING, {
        success: true,
        data: agentMessage,
        finishReason: 'stop',
        conversationUUID,
      });
      try {
        const routerAgents = this.routerAgents.get(conversationUUID);
        const searchAgent =
          routerAgents.functionToAgents[
            AgentFunctionByType[AgentType.SEQUENCE_SEARCH]
          ] || null;
        if (searchAgent) {
          await searchAgent.saveMessage(`${JSON.stringify(option)}`, Role.User);
          await searchAgent.saveMessage(
            `${JSON.stringify(res.data)}`,
            Role.Assistant,
          );
        }
      } catch (error) {
        this.logger.debug(error);
      }
    }
    return res;
  }

  async updateNcbiSearchMessage({
    messageId,
    option,
    optionRes,
  }: {
    messageId: number;
    selectedOptions?: Array<string>;
    option?: Record<any, any>;
    optionRes?: Record<any, any>;
  }) {
    const submitFromMessage = await this.dataSource.manager.findOne(
      MessageEntity,
      {
        where: {
          agentType: AgentType.SEQUENCE_SEARCH,
          id: messageId,
        },
      },
    );
    submitFromMessage.optionInfo = {
      ...submitFromMessage.optionInfo,
      ...option,
      optionRes,
    };
    this.logger.debug('-------->submitFromMessage.optionInfo.optionRes');
    this.logger.debug(submitFromMessage.optionInfo.optionRes);
    await this.dataSource.manager.update(
      MessageEntity,
      messageId,
      submitFromMessage,
    );
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
      await entityManager.delete(AgentMessageEntity, {
        conversationUUID,
      });

      await entityManager.delete(MessageEntity, {
        conversationUUID,
      });

      await entityManager.delete(ConversationEntity, {
        uuid: conversationUUID,
      });
      // todo:删除agents表
      // todo:删除devices表
      // todo:删除llms表
      // todo:删除tool表
    });
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
    // const chatHistory = history.map<ChatMessage>((item) => {
    //   return { role: item.role, content: item.content };
    // });
    // const messages = [
    //   ...this.initMessages,
    //   ..._.slice(chatHistory, -10), // TODO 应该计算token能覆盖到前面几条记录，而不是直接取最后10条
    //   { role: 'user', content: text },
    // ];
    // // this.logger.log(`messages: ${JSON.stringify(messages)}`);
    //
    // const functions = this._getFunctions();
    // // this.logger.log(`functions: ${JSON.stringify(functions)}`);
    //
    // const events = await this.aiClient.listChatCompletions(
    //   process.env.OPENAI_API_ENGINE,
    //   messages,
    //   {
    //     functionCall: 'auto',
    //     functions: functions,
    //   },
    // );
    //
    // await this._handleChatCompletions(
    //   client,
    //   text,
    //   assistantMessage,
    //   messages,
    //   events,
    // );
  }

  async updateUploadFileMessage(userMessage: MessageEntity, messageID: number) {
    await this.dataSource.manager.save(userMessage);
    await this.dataSource.manager.delete(MessageEntity, messageID);
  }
  setClient(client: Socket) {
    this.wsClient = client;
  }
  setOtDevices(otDevices: OtDevices) {
    this.otDevices = otDevices;
  }

  setNcbiSearch(ncbiSearch: NcbiSearch) {
    this.ncbiSearch = ncbiSearch;
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
