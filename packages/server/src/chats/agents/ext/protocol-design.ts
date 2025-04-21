import * as WebSocket from 'ws';
import { Operation } from './ncbi-search';
import { Logger } from '@nestjs/common';
import { DataSource } from 'typeorm';
import { MessageEntity } from 'src/chats/entities/message.entity';
import { AgentMessageEntity } from '../../entities/agent-message.entity';
import { CacheService } from 'src/cache/cache.service';
import { ChatMessage } from '@azure/openai';
import { AgentType } from '@xpcr/common';

// Protocol Design

// Protocol Design类型
export enum ProtocolDesignFunctionType {
  auto_protocol_design = 'auto_protocol_design',
  template_protocol_design = 'template_protocol_design',
}

// Protocol Design 入参
export interface ProtocolDesignInfo {
  stage: number;
  protocol_design_type: ProtocolDesignFunctionType;
  operations?: Operation[];
  data?: object;
}

// Protocol Design 出参
export interface ProtocolDesignResultInfo {
  responses: {
    response: string;
    stage: number;
    protocol_design_type: ProtocolDesignFunctionType;
    state: 'continue' | 'stop';
    operations?: Operation[];
    data?: object;
  };
}

export interface ProtocolDesignResult {
  content: ProtocolDesignResultInfo;
  state?: string;
}

export class ProtocolDesign {
  logger = new Logger(ProtocolDesign.name);
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
      `${AgentType.PROTOCOL_DESIGN}:${this.conversationUUID}`,
      false,
    );
    this.initConnect();
  }

  private initConnect() {
    this.ip = process.env.PROTOCOL_HOST;
    this.port = process.env.PROTOCOL_PORT;
    this.path = process.env.PROTOCOL_PATH;
    this.domain = process.env.PROTOCOL_DOMAIN;
    let url = '';
    if (this.domain) {
      url = `wss://${this.domain}${this.path}`;
    } else {
      url = `ws://${this.ip}:${this.port}${this.path}`;
    }
    this.logger.debug('protocol-design===>initConnect====>', url);
    this.ws = new WebSocket(url);

    this.ws.on('open', () => {
      this.logger.log(`connect success.`);
      this.connected = true;
    });

    this.ws.on('error', (error) => {
      this.logger.error(`error: ${error.toString()}`);
      this.sendResolve({
        responses: {
          response: `Protocol design error: ${error.toString()}`,
          state: 'stop',
          stage: 4,
        },
      } as ProtocolDesignResultInfo);
    });

    this.ws.on('close', () => {
      this.logger.log(`close`);
      this.cacheService.setObject(
        `${AgentType.PROTOCOL_DESIGN}:${this.conversationUUID}`,
        false,
      );
    });

    this.ws.on('message', (data) => {
      const text = data.toString();
      if (text == 'heartbeat') {
        this.logger.debug('Protocol Design heartbeat');
      } else {
        this.logger.log('== Protocol Design Result ==');
        this.logger.log(text);
        const protocolResult: ProtocolDesignResult = JSON.parse(text);
        this.cacheService.setObject(
          `${AgentType.PROTOCOL_DESIGN}:${this.conversationUUID}`,
          false,
        );
        if (typeof protocolResult.content == 'string') {
          const errorTip = `We apologize to the failure to help immediately;It seems that our assistants have encountered a little problem;Please try to check if there is an unsubmitted operation, our assistant will be happy to serve you.`;
          this.sendResolve(errorTip);
        } else {
          this.sendResolve(protocolResult.content);
        }
      }
    });
  }

  private sendResolve(content: ProtocolDesignResultInfo | string) {
    if (this.resolve) {
      this.resolve(content);
      this.resolve = null;
    }
  }

  async call(
    input: string,
    history: AgentMessageEntity[],
    optionInfo: ProtocolDesignInfo,
  ) {
    const generating = await this.cacheService.getObject(
      `${AgentType.PROTOCOL_DESIGN}:${this.conversationUUID}`,
    );
    if (generating) {
      return Promise.resolve<ProtocolDesignResultInfo>({
        responses: {
          response: `Protocol Design is running, please wait...`,
          stage: 2,
          protocol_design_type:
            ProtocolDesignFunctionType.template_protocol_design,
          state: 'continue',
        },
      });
    }

    this.cacheService.setObject(
      `${AgentType.PROTOCOL_DESIGN}:${this.conversationUUID}`,
      true,
    );

    const operationsObj = {};
    // 将operations里的用户提交的key解开放到外层
    if (optionInfo?.operations) {
      optionInfo?.operations.forEach((item) => {
        operationsObj[item.key] = item.value;
      });
    }

    const sendObj = {
      instruction: JSON.stringify({
        ...operationsObj,
        stage: optionInfo?.stage ?? 1,
        protocol_design_type:
          optionInfo?.protocol_design_type ??
          ProtocolDesignFunctionType.template_protocol_design,
        // 上一次的返回 不做修改
        data: optionInfo?.data ?? {},
      }),
      history,
      type: 0, // 无用字段
    };
    if (this.ws.readyState == WebSocket.CLOSED) {
      this.initConnect();
    }
    this.ws.send(JSON.stringify(sendObj));
    this.logger.log(
      `Protocol Design websocket send to backend: ${JSON.stringify(sendObj)}`,
    );

    return new Promise<ProtocolDesignResultInfo>((resolve) => {
      this.resolve = resolve;
    });
  }
}
