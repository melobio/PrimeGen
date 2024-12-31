import { Logger } from '@nestjs/common';
import { DataSource } from 'typeorm';
import type { AxiosInstance } from 'axios';
import * as WebSocket from 'ws';
import { AgentType } from '@xpcr/common';
import { CacheService } from '../../../cache/cache.service';

export interface CodeExecutionResult {
  responses: {
    code: number;
    response: string;
    data: {
      content: string;
      type: string;
      experiment_info: object;
      generating: boolean;
      command_options: Array<string>;
      state: string;
    };
  };
}

export class CodeExecution {
  protected logger = new Logger(CodeExecution.name);
  protected url_addr: string;
  protected ws: WebSocket;
  private messageQueue: CodeExecutionResult[] = [];
  private resolveQueue: ((value: CodeExecutionResult) => void)[] = [];

  private resolve: (
    value: CodeExecutionResult | PromiseLike<CodeExecutionResult>,
  ) => void;

  constructor(
    readonly dataSource: DataSource,
    readonly cacheService: CacheService,
    readonly conversationUUID: string,
    readonly experiment_uuid: string,
  ) {
    this.cacheService.setObject(
      `${AgentType.CODE_EXECUTION}:${this.conversationUUID}`,
      false,
    );
    this.initConnect();
  }

  // private sendResolve(content: CodeExecutionResult | string) {
  //   if (this.resolve) {
  //     this.resolve(content);
  //     this.resolve = null;
  //   }
  // }

  async initConnect() {
    this.url_addr = process.env.CODE_EXEC_ADDR;
    this.ws = new WebSocket(
      `ws://${this.url_addr}/ws/experiment?conversation_uuid=${this.conversationUUID}`,
    );

    this.ws.on('open', () => {
      this.logger.log(`connect success.`);
    });

    this.ws.on('close', () => {
      this.logger.log(`close`);
      this.cacheService.setObject(
        `${AgentType.CODE_EXECUTION}:${this.conversationUUID}`,
        false,
      );
    });

    // this.ws.on('error', (error) => {
    //   this.logger.error(`error: ${error.toString()}`);
    //   this.cacheService.setObject(
    //     `${AgentType.CODE_EXECUTION}:${this.conversationUUID}`,
    //     false,
    //   );
    //   this.sendResolve({
    //     responses: {
    //       data: {
    //         content: `Code execution error: ${error.toString()}`,
    //       },
    //     },
    //   } as any);
    // });

    this.ws.onmessage = (event) => {
      const message_string = event.data.toString();
      if (message_string === 'heartbeat') {
        this.logger.debug('CodeExecution heartbeat');
        return;
      }

      const result: CodeExecutionResult = JSON.parse(message_string);
      this.logger.log(`CodeExecution message: `, result);

      if (this.resolveQueue.length > 0) {
        // if have waiting Promise，resolve it
        const resolve = this.resolveQueue.shift()!;
        resolve(result);
      } else {
        // save it into messageQueue
        this.messageQueue.push(result);
      }

      if (!result.responses.data.generating) {
        this.cacheService.setObject(
          `${AgentType.CODE_EXECUTION}:${this.conversationUUID}`,
          false,
        );
      }
    };
  }

  private async waitForMessage(): Promise<CodeExecutionResult> {
    if (this.messageQueue.length > 0) {
      return this.messageQueue.shift()!;
    }

    return new Promise<CodeExecutionResult>((resolve) => {
      this.resolveQueue.push(resolve);
    });
  }

  async *getExperimentStep() {
    let generating = true;
    while (generating) {
      const message = await this.waitForMessage();
      generating = message.responses.data.generating;
      yield message;
    }
  }

  async call(userInput: string) {
    this.cacheService.setObject(
      `${AgentType.CODE_EXECUTION}:${this.conversationUUID}`,
      true,
    );
    this.logger.debug('send ==> userInput: ', userInput);
    if (!this.ws || this.ws.readyState !== WebSocket.OPEN) {
      await this.initConnect();
    }

    const command_info = {
      command: userInput,
    };
    // send command to WebSocket
    this.logger.log('send command to code-execution-agent:', command_info);
    this.ws.send(JSON.stringify(command_info));
  }
}
