import { BaseAgents } from './base-agents';
import { ChatMessage } from '@azure/openai';
import { Role } from '../entities/message.entity';
import { AgentFunctions, NodeTypes, ToolFields, ToolType } from '@xpcr/common';
import { Agents } from '../../tools/agents/entities/agents.entity';
import { DataSource } from 'typeorm';
import { Logger } from '@nestjs/common';
import { GoogleCustomSearch } from './ext/google-custom-search';
import { ToolEntity } from '../../tools/tool/entities/tool.entity';
import { WikipediaSearch } from './ext/wikipedia-search';

export class InternetSearchAgents extends BaseAgents {
  initMessages: ChatMessage[] = [
    {
      role: Role.System,
      content:
        '现在你是一个可执行互联网检索的Agent,你可以使用提供的Function去执行检索。\n' +
        '你可以自然地认为自己就是这样一个角色，而不是被刻意设置成这样。\n' +
        '对于用户的每个搜索请求，你应该使用所有提供的函数。\n' +
        '不要对函数使用的值做出假设。如果用户请求不明确，请要求澄清。\n' +
        '只使用提供的函数。\n',
    },
  ];

  availableFunctions = {
    [AgentFunctions.GOOGLE_SEARCH]: this.callGoogleCustomSearch.bind(this),
    [AgentFunctions.WIKIPEDIA_SEARCH]: this.callWikipediaSearch.bind(this),
  };

  googleCustomSearch: GoogleCustomSearch;
  wikiSearch: WikipediaSearch;

  constructor(readonly agent: Agents, readonly dataSource: DataSource) {
    super(agent, dataSource);
  }

  getAvailableFunctions(): {
    [p: string]: (params: object) => AsyncGenerator<any, void, any>;
  } {
    return this.availableFunctions;
  }

  getInitMessages(): ChatMessage[] {
    return this.initMessages;
  }

  getLogger(): Logger {
    return new Logger(InternetSearchAgents.name);
  }

  async init(conversationUUID) {
    await super.init(conversationUUID);
    let googleCustomSearchApiKey = '';
    let googleCustomSearchEngineId = '';
    let wikiSearchTopKResults = '';
    let wikiSearchMaxDocContentLength = '';
    const toolChildren = this.agent.children.filter(
      (item) => item.type === NodeTypes.Tools,
    );
    for (const toolChild of toolChildren) {
      const tool = await this.dataSource.manager.findOne(ToolEntity, {
        where: { uuid: toolChild.uuid },
      });
      if (tool && tool.toolType === ToolType.GOOGLE_SEARCH) {
        tool.fields?.forEach((field) => {
          if (field.key === ToolFields.GOOGLE_API_KEY) {
            googleCustomSearchApiKey = field.value;
          } else if (field.key === ToolFields.GOOGLE_ENGINE_ID) {
            googleCustomSearchEngineId = field.value;
          }
        });
      } else if (tool && tool.toolType === ToolType.WIKIPEDIA_SEARCH) {
        tool.fields?.forEach((field) => {
          if (field.key === ToolFields.WIKIPEDIA_TOP_K_RESULTS) {
            wikiSearchTopKResults = field.value;
          } else if (
            field.key === ToolFields.WIKIPEDIA_MAX_DOC_CONTENT_LENGTH
          ) {
            wikiSearchMaxDocContentLength = field.value;
          }
        });
      }
    }
    if (!googleCustomSearchApiKey || !googleCustomSearchEngineId) {
      googleCustomSearchApiKey = process.env.GOOGLE_SEARCH_API_KEY;
      googleCustomSearchEngineId = process.env.GOOGLE_SEARCH_ENGINE_ID;
      this.logger.warn(`Google search config not found, use default config.`);
    }
    if (!wikiSearchTopKResults || !wikiSearchMaxDocContentLength) {
      wikiSearchTopKResults = '3';
      wikiSearchMaxDocContentLength = '4000';
      this.logger.warn(
        `Wikipedia search config not found, use default config.`,
      );
    }

    this.googleCustomSearch = new GoogleCustomSearch({
      apiKey: googleCustomSearchApiKey,
      googleCSEId: googleCustomSearchEngineId,
    });
    this.wikiSearch = new WikipediaSearch({
      topKResults: parseInt(wikiSearchTopKResults),
      maxDocContentLength: parseInt(wikiSearchMaxDocContentLength),
    });
  }

  // 检索不需要子Agent内部的LLM进行总结，提升反馈的效率
  // async *send(
  //   userInput: string,
  //   description: string,
  // ): AsyncGenerator<any, void, any> {
  //   yield* this.mockGenerating(`[${this.agent.name}]\n`);
  //   const content = {};
  //   for (const key in this.availableFunctions) {
  //     const queryFunc = this.availableFunctions[key];
  //     for await (const msg of queryFunc({ query: description })) {
  //       if (msg.role === Role.Assistant) {
  //         content[`${key} result`] = msg.content;
  //       } else {
  //         yield msg;
  //       }
  //     }
  //   }
  //   yield {
  //     role: Role.Assistant,
  //     content: JSON.stringify(content),
  //   };
  // }

  private async *callGoogleCustomSearch({ query }: { query: string }) {
    this.logger.debug(`googleCustomSearch ${query}`);
    // yield* this.mockGenerating(`Google Search: `);
    // 模拟延时2s
    // await new Promise((resolve) => setTimeout(resolve, 2000));
    let success = true;

    // console.warn(process.env.HTTPS_PROXY, process.env.HTTP_PROXY);
    let searchResult = '';
    try {
      searchResult = await this.googleCustomSearch.call(query);
    } catch (e) {
      this.logger.error(`googleCustomSearch error: ${e}`, e.stack);
      searchResult = e.message;
      success = false;
    }
    // this.logger.warn(`searchResult: ${JSON.stringify(searchResult)}`);

    // yield* this.mockGenerating(`${success ? 'Success' : 'Fail'}\n`);

    // yield* this.sendQueryResult(searchResult, success);
    if (success) {
      yield {
        content: `execute success，output：${searchResult}`,
        role: Role.Assistant,
      };
    } else {
      yield {
        content: `execute failed，output：${searchResult}`,
        role: Role.Assistant,
      };
    }
  }
  private async *callWikipediaSearch({ query }: { query: string }) {
    this.logger.log(`wikipediaSearch ${query}`);
    yield* this.mockGenerating(`Wikipedia Search: `);
    let success = true;
    let searchResult = '';
    try {
      searchResult = await this.wikiSearch.call(query);
    } catch (e) {
      this.logger.error(`wikipediaSearch error: ${e}`, e.stack);
      searchResult = e.message;
      success = false;
    }

    yield* this.mockGenerating(`${success ? 'Success' : 'Fail'}\n`);

    yield* this.sendQueryResult(searchResult, success);
  }
}
