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

        // delete agents table
        await entityManager.delete(Agents, {
          conversationUUID,
        });
        // delete devices table
        await entityManager.delete(DevicesEntity, {
          conversationUUID,
        });
        // delete llms table
        await entityManager.delete(LlmsEntity, {
          conversationUUID,
        });
        // delete tool table
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
    //
  }

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

  // call Agent
  private async *_callFunction(functionCall: FunctionCall) {
    //
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
    // all history messages
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
}
