import { Logger } from '@nestjs/common';
import * as WebSocket from 'ws';
import { NCBIResultContent } from './ncbi-search';
import { AgentMessageEntity } from '../../entities/agent-message.entity';
import { DataSource } from 'typeorm';
import { MessageEntity } from 'src/chats/entities/message.entity';
import { AgentType } from '@xpcr/common';
import { ChatMessage } from '@azure/openai';

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

  protected generating = false;
  protected resolve?: (value: any) => void;
  protected ws: WebSocket;

  protected connected = false;
  constructor(
    readonly dataSource: DataSource,
    readonly conversationUUID: string,
  ) {
    this.ip = process.env.PRIMER_HOST;
    this.port = process.env.PRIMER_PORT;
    this.path = process.env.PRIMER_PATH;

    const url = `ws://${this.ip}:${this.port}${this.path}`;
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
      this.generating = false;
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
        this.generating = false;
        if (primerResult.state === 'stop') {
          // this.resolve(this.current.content);
          this.sendResolve(primerResult.content);
        }
      }
    });
  }

  private sendResolve(content: PrimerResultContent) {
    if (this.resolve) {
      this.logger.log(`Primer Design result: ${content}`);
      this.resolve(content);
      this.resolve = null;
    }
  }

  async call(
    input: string,
    history: AgentMessageEntity[],
    lastPrimerMsg: Record<string, any> = {},
  ) {
    if (this.generating) {
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
    this.generating = true;
    const chatHistory = await this.dataSource.manager.find(MessageEntity, {
      where: { conversationUUID: this.conversationUUID },
      order: { createTime: 'asc' },
    });
    const chatHistoryArr = chatHistory.map<ChatMessage>((item) => {
      return { role: item.role, content: item.content };
    });
    const sendObj = {
      instruction: JSON.stringify({
        stage: 1,
        ...lastPrimerMsg.optionInfo,
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
