import { Logger } from '@nestjs/common';
import { DataSource } from 'typeorm';
import { ChatMessage } from '@azure/openai';
import type { AxiosInstance } from 'axios';
import * as WebSocket from 'ws';
import * as fs from 'fs';
import * as path from 'path';
import axios from 'axios';
import { AgentType, AgentFunctions } from '@xpcr/common';

import { BaseAgents } from './base-agents';
import { ConversationService } from '../conversation.service';
import { CacheService } from '../../cache/cache.service';
import { Agents } from '../../tools/agents/entities/agents.entity';
import { Role } from '../entities/message.entity';
import { assert } from 'console';
import { CodeExecution, CodeExecutionResult } from './ext/code-execution';

export class CodeExecutionAgent extends BaseAgents {
  // attributes
  initMessages: ChatMessage[] = [];
  codeExecution: CodeExecution;
  // LLM available functions for function calling
  availableFunctions = {
    [AgentFunctions.CODE_EXECUTION]: this.send.bind(this),
  };
  url_addr: string;
  http_client: AxiosInstance;
  experiment_uuid: string;

  protected resolve?: (value: any) => void;

  // Initialize the CodeExecutionAgent
  constructor(
    readonly agent: Agents,
    readonly dataSource: DataSource,
    readonly cacheService: CacheService,
    private readonly chatsService: ConversationService,
  ) {
    super(agent, dataSource);
    this.url_addr = process.env.CODE_EXEC_ADDR;
    this.http_client = axios.create({
      baseURL: `http://${this.url_addr}`,
      timeout: 10000,
    });
  }

  // Override Functions
  getAvailableFunctions(): {
    [p: string]: (params: object) => AsyncGenerator<any, void, any>;
  } {
    return this.availableFunctions;
  }

  getInitMessages(): ChatMessage[] {
    return this.initMessages;
  }

  getLogger(): Logger {
    return new Logger(CodeExecutionAgent.name);
  }
  // End Override Functions

  async get_protocol_design_result() {
    assert(this.chatsService.protocolDesigns.has(this.conversationUUID));
    // this chat history message, filter only protocol design message
    const history = await this.chatsService.findAllMessages(
      this.conversationUUID,
    );
    // protocol design result
    const protocolDesignMessages = history.filter(
      (message) =>
        message.role === Role.Assistant &&
        message.agentType === AgentType.PROTOCOL_DESIGN &&
        message.conversationUUID === this.conversationUUID &&
        message.optionInfo &&
        message.optionInfo?.state == 'stop',
    );
    if (protocolDesignMessages.length < 1) {
      return {};
    }
    const lastOptionInfo = protocolDesignMessages[0].optionInfo;
    const protocolDesignResult = {
      conversation_uuid: this.conversationUUID,
      protocol_path: lastOptionInfo.data?.json_file,
      panel_list: lastOptionInfo.data?.layout_info,
    };
    return protocolDesignResult;
  }

  async init(conversationUUID) {
    const codeExecutionAgent = fs.readFileSync(
      path.join(process.cwd(), 'assets', 'prompts', 'code-execution-agent.txt'),
    );
    this.initMessages = [
      {
        role: Role.System,
        content: String(codeExecutionAgent),
      },
    ];
    await super.init(conversationUUID);
  }

  // Send data to python backend /ws/experiment websocket api
  async *send({ userInput, description }) {
    // get experiment id
    if (!this.experiment_uuid) {
      this.logger.error('experiment_uuid is null, try to get experiment_uuid');
      const protocolDesignResult = await this.get_protocol_design_result();
      if (protocolDesignResult) {
        const resp = await this.http_client.post(
          '/api/conversation/experiment',
          protocolDesignResult,
        );
        this.experiment_uuid = resp.data.experiment_uuid;
      }
    }
    // object to communication with python
    this.codeExecution = new CodeExecution(
      this.dataSource,
      this.cacheService,
      this.conversationUUID,
      this.experiment_uuid,
    );

    this.logger.debug(
      `send ==> userInput: ${userInput}, description: ${description}`,
    );
    await this.codeExecution.call(description);
    console.log(`Send to code-execution-agent: ${description}`);

    const msgs = [];
    for await (const resp of this.codeExecution.getExperimentStep()) {
      if (resp && resp.responses) {
        this.logger.log(`resp.responses:`, resp.responses);
        if (msgs.length > 3) {
          msgs.shift();
        }
        this.logger.log('CurrentMessagesLength:', msgs.length);
        msgs.push(resp.responses.data.content);
        const content = msgs.join('\n');
        await this.saveMessage(`${userInput} ${description}`, Role.User);
        await this.saveMessage(content, Role.Assistant);
        yield {
          role: Role.Assistant,
          content: content,
          optionInfo: resp.responses,
        };
      }
    }
  }
}
