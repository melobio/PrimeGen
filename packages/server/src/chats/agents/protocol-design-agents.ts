import { BaseAgents } from './base-agents';
import { ChatMessage } from '@azure/openai/types/src';
import { AgentFunctions, AgentType } from '@xpcr/common';
import { Logger } from '@nestjs/common';
import { Agents } from '../../tools/agents/entities/agents.entity';
import { DataSource } from 'typeorm';
import { ConversationService } from '../conversation.service';
import { CacheService } from 'src/cache/cache.service';
import * as fs from 'fs';
import * as path from 'path';
import { MessageEntity, Role } from '../entities/message.entity';
import {
  ProtocolDesignInfo,
  ProtocolDesign,
  ProtocolDesignResultInfo,
} from './ext/protocol-design';

const createFileListString = (key, responses) => {
  const fileList = responses[key] || [];
  if (fileList.filter((a) => a).length === 0) return '';
  const keyName = key?.replaceAll('_', ' ').toUpperCase() || '';
  const fileListLinks = fileList.map((p) => {
    const url = path.join('/v1/files/xsearchdev', p);
    return `[${path.basename(url)}](${url.replaceAll(' ', '%20')})`;
  });

  return `**${keyName}**:\n${fileListLinks.join('\n')};\n`;
};

export class ProtocolDesignAgents extends BaseAgents {
  initMessages: ChatMessage[] = [];
  protocolDesign: ProtocolDesign;
  availableFunctions = {
    [AgentFunctions.DESIGN_PROTOCOL]: this.callProtocolDesign.bind(this),
  };

  constructor(
    readonly agent: Agents,
    readonly dataSource: DataSource,
    readonly cacheService: CacheService,
    readonly chatsService: ConversationService,
) {
    super(agent, dataSource);
    this.needSummarize = false;
  }

  async init(conversationUUID) {
    const protocolDesignAgent = fs.readFileSync(
      path.join(
        process.cwd(),
        'assets',
        'prompts',
        'protocol-design-agent.txt',
      ),
    );
    this.initMessages = [
      {
        role: Role.System,
        content: String(protocolDesignAgent),
      },
    ];
    await super.init(conversationUUID);
    this.protocolDesign = new ProtocolDesign(
      this.dataSource,
      this.cacheService,
      this.conversationUUID,
    );
    this.chatsService.setProtocolDesign(
      this.protocolDesign,
      this.conversationUUID,
    );
  }


  async *send({
    userInput,
    description,
    optionInfo,
  }: {
    userInput: string;
    description?: string;
    optionInfo?: ProtocolDesignInfo;
  }): AsyncGenerator<any, void, any> {
    // Send User Input to Protocol Design Agent
    this.logger.debug(
      `send ==> userInput: ${userInput} description: ${description} optionInfo: ${JSON.stringify(
        optionInfo,
      )}`,
    );

    let content = '';
    let newOptionInfo = {
      state: 'stop',
      // stage: 1,
      operations: [],
    };
    for (const key in this.availableFunctions) {
      const queryFunc = this.availableFunctions[key];
      for await (const msg of queryFunc({
        query: `${userInput} ${description}`,
        optionInfo,
      })) {
        if (msg?.optionInfo) {
          newOptionInfo = msg?.optionInfo;
        }
        if (msg.role === Role.Assistant) {
          content += msg.content;
        } else {
          yield msg;
        }
      }
    }
    yield {
      role: Role.Assistant,
      content: content,
      optionInfo: newOptionInfo,
    };
    await this.saveMessage(`${userInput} ${description}`, Role.User);
    await this.saveMessage(content, Role.Assistant);
  }

  private async *callProtocolDesign({
    query,
    optionInfo,
  }: {
    query: string;
    optionInfo?: ProtocolDesignInfo;
  }) {
    let success = true;
    let protocolResultString = '';
    let protocolResultContentData;
    let protocolResultContentResponse: ProtocolDesignResultInfo;
    let protocolResultContent;

    try {
      protocolResultContentResponse = await this.protocolDesign.call(
        query,
        this.history,
        optionInfo,
      );

      if (typeof protocolResultContent == 'string') {
        protocolResultString = `We apologize to the failure to help immediately;It seems that our assistants have encountered a little problem;Please try to check if there is an unsubmitted operation, our assistant will be happy to serve you.`;
      } else {
        protocolResultContent =
          protocolResultContentResponse?.responses || null;
        protocolResultContentData = protocolResultContent.data;
        protocolResultString = protocolResultContent?.response || '';

        if (protocolResultContentData?.new_code) {
          protocolResultString += '\n';
          protocolResultString += createFileListString(
            'new_code',
            protocolResultContentData,
          );
        }

        if (protocolResultContentData?.json_file) {
          protocolResultString += '\n';
          protocolResultString += createFileListString(
            'json_file',
            protocolResultContentData,
          );
        }

        if (protocolResultContentData?.layout_info) {
          protocolResultString += '\n';
          protocolResultString += createFileListString(
            'layout_info',
            protocolResultContentData,
          );
        }
      }
    } catch (e) {
      this.logger.error(`designProtocol error: ${e}`, e.stack);
      protocolResultString = e.message;
      success = false;
    }

    yield {
      optionInfo: protocolResultContent,
      role: Role.Assistant,
      content: protocolResultString,
    };
  }

  getLogger(): Logger {
    return new Logger(ProtocolDesignAgents.name);
  }
  getInitMessages(): ChatMessage[] {
    return this.initMessages;
  }
  getAvailableFunctions(): {
    [p: string]: (params: object) => AsyncGenerator<any, void, any>;
  } {
    return this.availableFunctions;
  }
}
