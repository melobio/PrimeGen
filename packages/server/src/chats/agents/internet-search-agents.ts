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
import { CacheService } from 'src/cache/cache.service';

export class InternetSearchAgents extends BaseAgents {
  initMessages: ChatMessage[] = [
    {
      role: Role.System,
      content:
        'Now that you are an Agent that can perform Internet searches, you can use the provided Function to perform the searches. \n' +
        'You can naturally assume that you are such a role, rather than being deliberately set up as such. \n' +
        'For each search request by a user, you should use all provided Functions. \n' +
        "Don't make assumptions about the values used by the functions. If a user request is unclear, ask for clarification. \n" +
        'Use only the provided functions. \n',
    },
  ];

  availableFunctions = {
    [AgentFunctions.GOOGLE_SEARCH]: this.callGoogleCustomSearch.bind(this),
    [AgentFunctions.WIKIPEDIA_SEARCH]: this.callWikipediaSearch.bind(this),
  };

  googleCustomSearch: GoogleCustomSearch;
  wikiSearch: WikipediaSearch;

  constructor(
    readonly agent: Agents,
    readonly dataSource: DataSource,
    readonly cacheService: CacheService,
  ) {
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

  private async *callGoogleCustomSearch({ query }: { query: string }) {
    this.logger.debug(`googleCustomSearch ${query}`);
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
