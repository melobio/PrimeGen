import { Logger } from '@nestjs/common';
import * as WebSocket from 'ws';
import { NCBIResultContent, SnpPrimerDesignInfo } from './ncbi-search';
import { AgentMessageEntity } from '../../entities/agent-message.entity';
import { DataSource } from 'typeorm';
import { MessageEntity } from 'src/chats/entities/message.entity';
import { AgentType } from '@xpcr/common';
import { ChatMessage } from '@azure/openai';
import { CacheService } from 'src/cache/cache.service';

export interface PrimerResultContent {
  responses: {
    response: string;
    operations: any[];
    state: 'stop' | 'continue';
    primer_type: PrimerType;
    stage: number;
    primer_design_prompt?: string;
    primer_design_dict?: Record<string, any>;
    primer?: string[];
  };
  history_conversation: Array<Record<string, string>>;
}

type PrimerType = 'snp_primer_design_type';
export interface PrimerResult {
  content: PrimerResultContent;
  type: number;
  state?: string;
}

export class PrimerDesign {
  logger = new Logger(PrimerDesign.name);
  protected ip: string;
  protected port: string;
  protected path: string;
  protected domain: string;
  protected resolve?: (value: any) => void;
  protected ws: WebSocket;

  protected connected = false;
  constructor(
    readonly dataSource: DataSource,
    readonly cacheService: CacheService,
    readonly conversationUUID: string,
  ) {
    this.cacheService.setObject(
      `${AgentType.PRIMER_DESIGN}:${this.conversationUUID}`,
      false,
    );
    this.ip = process.env.PRIMER_HOST;
    this.port = process.env.PRIMER_PORT;
    this.path = process.env.PRIMER_PATH;
    this.domain = process.env.PRIMER_DOMAIN;
    let url = '';
    if (this.domain) {
      url = `wss://${this.domain}${this.path}`;
    } else {
      url = `ws://${this.ip}:${this.port}${this.path}`;
    }
    this.logger.debug('primer-design===>initConnect====>', url);
    this.ws = new WebSocket(url);

    this.ws.on('open', () => {
      this.logger.log(`connect success.`);
      this.connected = true;
    });

    this.ws.on('error', (error) => {
      this.logger.error(`error: ${error.toString()}`);
      this.sendResolve({
        responses: {
          response: `Primer design error: ${error.toString()}`,
        },
      } as PrimerResultContent);
    });

    this.ws.on('close', () => {
      this.logger.log(`close`);
      this.cacheService.setObject(
        `${AgentType.PRIMER_DESIGN}:${this.conversationUUID}`,
        false,
      );
    });

    this.ws.on('message', (data) => {
      const text = data.toString();
      // this.logger.log(`message: ${text}`);
      if (text === 'heartbeat') {
        this.logger.log(`primer design heartbeat`);
      } else {
        this.logger.debug('== Primer Design Result ==');
        this.logger.log(text);
        const primerResult: PrimerResult = JSON.parse(text);
        this.cacheService.setObject(
          `${AgentType.PRIMER_DESIGN}:${this.conversationUUID}`,
          false,
        );
        if (primerResult.state === 'stop') {
          this.sendResolve(primerResult.content);
        }
      }
    });
  }

  private sendResolve(content: PrimerResultContent) {
    if (this.resolve) {
      this.logger.log(`Primer Design result: ${JSON.stringify(content)}`);
      this.resolve(content);
      this.resolve = null;
    }
  }

  async call(
    input: string,
    history: AgentMessageEntity[],
    lastPrimerMsg: Record<string, any> = {},
    optionInfo: SnpPrimerDesignInfo,
  ) {
    const generating = await this.cacheService.getObject(
      `${AgentType.PRIMER_DESIGN}:${this.conversationUUID}`,
    );
    if (generating) {
      return Promise.resolve<PrimerResultContent>({
        responses: {
          response: `Primer Design is running, please wait...`,
        },
      } as PrimerResultContent);
    }
    const searchAgentMessages = await this.dataSource.manager.find(
      MessageEntity,
      {
        where: {
          agentType: AgentType.SEQUENCE_SEARCH,
          conversationUUID: this.conversationUUID,
        },
      },
    );
    const searchLastoptionInfo = searchAgentMessages.find(
      (item) => item.optionInfo && item.optionInfo.state == 'stop',
    ).optionInfo;
    const search_responses = searchLastoptionInfo || {};
    this.cacheService.setObject(
      `${AgentType.PRIMER_DESIGN}:${this.conversationUUID}`,
      true,
    );
    const chatHistory = await this.dataSource.manager.find(MessageEntity, {
      where: { conversationUUID: this.conversationUUID },
      order: { createTime: 'asc' },
    });
    const chatHistoryArr = chatHistory.map<ChatMessage>((item) => {
      return { role: item.role, content: item.content };
    });
    const operationsObj = {};
    if (optionInfo?.operations) {
      optionInfo?.operations.forEach((item) => {
        operationsObj[item.key] = item.value;
      });
    }
    const sendObj = {
      instruction: JSON.stringify({
        stage: 1,
        ...lastPrimerMsg.optionInfo,
        ...operationsObj,
        conversation: chatHistoryArr,
        search_responses,
      }),
      history,
      type: 0, // 无用字段
    };
    this.ws.send(JSON.stringify(sendObj));
    this.logger.log(`send: ${JSON.stringify(sendObj)}`);

    return new Promise<PrimerResultContent>((resolve) => {
      this.resolve = resolve;
    });
  }
}
