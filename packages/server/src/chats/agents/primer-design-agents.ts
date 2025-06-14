import { BaseAgents } from './base-agents';
import { ChatMessage } from '@azure/openai/types/src';
import { AgentFunctions, AgentType } from '@xpcr/common';
import { Logger } from '@nestjs/common';
import { Agents } from '../../tools/agents/entities/agents.entity';
import { DataSource } from 'typeorm';
import * as fs from 'fs';
import * as path from 'path';
import { MessageEntity, Role } from '../entities/message.entity';
import { PrimerDesign, PrimerResultContent } from './ext/primer-design';
import { CacheService } from 'src/cache/cache.service';
import { ConversationService } from '../conversation.service';
import { SnpPrimerDesignInfo } from './ext/ncbi-search';
const mockResult = `
  根据你的输入"{{INPUT}}",
  设计出的引物如下:
  
  primer1.csv
  primer2.csv
  primer3.csv
  primer4.csv
`;

export class PrimerDesignAgents extends BaseAgents {
  initMessages: ChatMessage[] = [];
  primerDesign: PrimerDesign;
  availableFunctions = {
    [AgentFunctions.DESIGN_PRIMER]: this.designPrimer.bind(this),
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
    const primerDesignAgent = fs.readFileSync(
      path.join(process.cwd(), 'assets', 'prompts', 'primer-design-agent.txt'),
    );
    this.initMessages = [
      {
        role: Role.System,
        content: String(primerDesignAgent),
      },
    ];
    await super.init(conversationUUID);
    this.primerDesign = new PrimerDesign(
      this.dataSource,
      this.cacheService,
      this.conversationUUID,
    );
    this.chatsService.setPrimerDesign(this.primerDesign, this.conversationUUID);
  }

  // 引物设计不需要子Agent内部的LLM进行总结，提升反馈的效率
  async *send({
    userInput,
    description,
    optionInfo,
  }: {
    userInput: string;
    description?: string;
    optionInfo?: SnpPrimerDesignInfo;
  }): AsyncGenerator<any, void, any> {
    // yield* this.mockGenerating(`[${this.agent.name}]\n`);
    let content = '';
    let newOptionInfo = { state: 'stop' };
    for (const key in this.availableFunctions) {
      const queryFunc = this.availableFunctions[key];
      for await (const msg of queryFunc({
        query: `${userInput} ${description}`,
        optionInfo,
      })) {
        if (msg?.optionInfo) {
          newOptionInfo = msg?.optionInfo;
        }
        // 将用户完整的输入传给Seq Agent
        if (msg.role === Role.Assistant) {
          content += msg.content;
        } else {
          yield msg;
        }
      }
    }
    yield {
      optionInfo: newOptionInfo,
      role: Role.Assistant,
      content: content,
    };
    await this.saveMessage(`${userInput} ${description}`, Role.User);
    await this.saveMessage(content, Role.Assistant);
  }

  private async *designPrimer({
    query,
    optionInfo,
  }: {
    query: string;
    optionInfo: SnpPrimerDesignInfo;
  }) {
    this.logger.debug(`designPrimer ${query}`);
    let success = true;
    let primerResult = '';
    let primerResultContent: PrimerResultContent;
    let resultOptionInfo: SnpPrimerDesignInfo;
    let primerResultContentResponses;
    if (!optionInfo) {
      optionInfo = {
        operations: [],
        state: 'continue',
        primer_type: '',
        stage: 2,
        primer_design_prompt: '',
        primer_design_dict: {},
      };
    }
    try {
      const allMessages = await this.dataSource.manager.find(MessageEntity, {
        where: { conversationUUID: this.conversationUUID },
      });
      const primerMsgs: MessageEntity[] = allMessages.filter(
        (message) => message.agentType === AgentType.PRIMER_DESIGN,
      );
      const lastPrimerMsg: MessageEntity = primerMsgs.at(-2);
      primerResultContent = await this.primerDesign.call(
        query,
        this.history,
        lastPrimerMsg,
        optionInfo,
      );
      if (typeof primerResultContent == 'string') {
        const errorTip = `We apologize to the failure to help immediately;It seems that our assistants have encountered a little problem;Please try to check if there is an unsubmitted operation, our assistant will be happy to serve you.`;
        primerResult = errorTip;
      } else {
        primerResultContentResponses = primerResultContent?.responses || null;

        if (primerResultContentResponses) {
          resultOptionInfo = {
            operations: primerResultContentResponses?.operations || [],
            state: primerResultContentResponses?.state || '',
            primer_type: primerResultContentResponses?.primer_type || '',
            stage: primerResultContentResponses?.stage || 2,
            primer_design_prompt:
              primerResultContentResponses?.primer_design_prompt || '',
            primer_design_dict:
              primerResultContentResponses?.primer_design_dict || {},
          };
          primerResult = `${primerResultContentResponses.response}`;
        }
        const primerList = primerResultContentResponses?.data?.primer || [];
        const primerUrlList = primerList.map((p) => {
          const url = path.join('/v1/files/xsearchdev', p);
          return `[${path.basename(url)}](${url.replaceAll(' ', '%20')})`;
        });
        primerResult = `${
          primerResultContentResponses.response || ''
        }\n${primerUrlList.join('\n')}\n`;
      }
    } catch (e) {
      primerResult = e.message;
      success = false;
    }

    yield {
      optionInfo: resultOptionInfo,
      content: primerResult,
      role: Role.Assistant,
    };
  }

  getLogger(): Logger {
    return new Logger(PrimerDesignAgents.name);
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
