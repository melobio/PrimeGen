import { BaseAgents } from './base-agents';
import { ChatMessage } from '@azure/openai';
import { MiddleMessageType, Role } from '../entities/message.entity';
import { AgentFunctions } from '@xpcr/common';
import {
  NCBIResultContent,
  NcbiSearch,
  TargetGeneInfo,
  NCBIFunctionType,
} from './ext/ncbi-search';
import { Agents } from '../../tools/agents/entities/agents.entity';
import { DataSource } from 'typeorm';
import { Logger } from '@nestjs/common';
import * as fs from 'fs';
import * as path from 'path';
import { ConversationService } from '../conversation.service';

const mockResult = `
  根据你的输入"{{INPUT}}",
  我在NCBI上检索到了如下结果：
  
  [Candidatus Mycobacterium methanotrophicum_cds_from_genomic.fna](https://www.ncbi.nlm.nih.gov/nuccore/NC_017939.1?report=fasta&from=1&to=100)
  [Candidatus Mycobacterium wuenschmannii cds from_genomic.fna](https://www.ncbi.nlm.nih.gov/nuccore/NC_017939.1?report=fasta&from=1&to=100)
  [Mycobacterium abscessus cds from_ genomic.fna](https://www.ncbi.nlm.nih.gov/nuccore/NC_017939.1?report=fasta&from=1&to=100)
  [Mycobacterium_adipatum_cds_from_genomic.fna](https://www.ncbi.nlm.nih.gov/nuccore/NC_017939.1?report=fasta&from=1&to=100)
  [Mycobacterium agri cds from_genomic.fna](https://www.ncbi.nlm.nih.gov/nuccore/NC_017939.1?report=fasta&from=1&to=100)
  [Mycobacterium aichiense cds from _ genomic.fna](https://www.ncbi.nlm.nih.gov/nuccore/NC_017939.1?report=fasta&from=1&to=100)
  [Mycobacterium_algericum_cds_from_genomic.fna](https://www.ncbi.nlm.nih.gov/nuccore/NC_017939.1?report=fasta&from=1&to=100)
  [Mycobacterium alsense_ cds from_genomic.fna](https://www.ncbi.nlm.nih.gov/nuccore/NC_017939.1?report=fasta&from=1&to=100)
  [Mycobacterium alsiense _ cds from_ genomic.fna](https://www.ncbi.nlm.nih.gov/nuccore/NC_017939.1?report=fasta&from=1&to=100)
  [Mycobacterium alvei cds from genomic.fna](https://www.ncbi.nlm.nih.gov/nuccore/NC_017939.1?report=fasta&from=1&to=100)
  
  请选择一个或者多个结果
`;

export class SequenceSearchAgents extends BaseAgents {
  initMessages: ChatMessage[] = [];
  availableFunctions = {
    [AgentFunctions.NCBI_SEARCH]: this.callNCBISearch.bind(this),
  };

  ncbiSearch: NcbiSearch;

  constructor(
    readonly agent: Agents,
    readonly dataSource: DataSource,
    private readonly chatsService: ConversationService,
  ) {
    super(agent, dataSource);
    this.needSummarize = false;
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
    return new Logger(SequenceSearchAgents.name);
  }

  async init(conversationUUID) {
    const seqSearchAgent = fs.readFileSync(
      path.join(
        process.cwd(),
        'assets',
        'prompts',
        'sequence-search-agent.txt',
      ),
    );
    this.initMessages = [
      {
        role: Role.System,
        content: String(seqSearchAgent),
      },
    ];
    await super.init(conversationUUID);
    this.ncbiSearch = new NcbiSearch();
    this.chatsService.setNcbiSearch(this.ncbiSearch);
  }

  // 检索不需要子Agent内部的LLM进行总结，提升反馈的效率
  async *send(
    userInput: string,
    description: string,
  ): AsyncGenerator<any, void, any> {
    // yield* this.mockGenerating(`[${this.agent.name}]\n`);
    let content = '';

    let optionInfo = {
      options: [],
      cancer_list: [],
      selected_options: [],
      stage: 2,
      state: '',
      search_type: '',
    };
    let operation = [];
    for (const key in this.availableFunctions) {
      const conversationMessages = await this.chatsService.findAllMessages(
        this.conversationUUID,
      );
      const chatHistory = conversationMessages.map<ChatMessage>((item) => {
        return { role: item.role, content: item.content };
      });
      const queryFunc = this.availableFunctions[key];
      for await (const msg of queryFunc({
        query: `${userInput}`,
        conversation: chatHistory,
      })) {
        this.logger.debug(`callNCBISearch==userInput===>`);
        this.logger.debug(userInput);
        if (msg?.optionInfo) {
          optionInfo = msg?.optionInfo;
        }
        if (msg?.operation) {
          operation = msg?.operation;
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
      role: Role.Assistant,
      content: content,
      optionInfo,
      operation,
    };
    await this.saveMessage(`${userInput} ${description}`, Role.User);
    await this.saveMessage(content, Role.Assistant);
  }

  private async *callNCBISearch({ query }: { query: string }) {
    this.logger.debug(`callNCBISearch==query===> ${query}`);
    // todo:拦截 用户选项的答复 推断后的函数，然后拼接参数后再给agent发消息
    // yield* this.mockGenerating(`Sequence Search: `);
    let success = true;
    let searchResult = '';
    let options: Array<Record<string, any>> = [];
    let cancer_list: Record<string, any> = {};
    let data_list: Array<Record<string, any>> = [];
    let stage = null;
    let state = '';
    let search_type = '';
    let response = '';
    let upload_file = false;
    let species_identification_dict = {};
    let operations = [];
    try {
      const res: NCBIResultContent = await this.ncbiSearch.call(
        query,
        this.history,
      );
      const responses = res.responses;
      const target_flie_path = responses?.target_flie_path || [];
      stage = responses.stage;
      state = responses.state;
      response = responses.response;
      search_type = responses.search_type;
      if (search_type == NCBIFunctionType.cancer_type) {
        cancer_list = responses?.cancer_list || {};
        searchResult = `${response}: \n`;
      } else if (
        search_type == NCBIFunctionType.genetic_disorder_type ||
        search_type == NCBIFunctionType.protein_mutation_type
      ) {
        // 遗传疾病/蛋白质突变
        if (responses?.options && responses.options.length > 0) {
          data_list = responses?.data_list || [];
          if (Array.isArray(responses.options)) {
            options = responses.options?.map((item, index) => ({
              ...item,
              index,
            }));
            searchResult = `${response}: \n`;
          } else {
            searchResult = `${response}: \n ${options}`;
          }
        }
      } else if (
        // 物种鉴定
        responses.search_type == NCBIFunctionType.species_identification_type
      ) {
        this.logger.debug('物种鉴定=================>');
        upload_file = responses.upload_file;
        operations = responses.operations;
        species_identification_dict = responses?.species_identification_dict;
        if (Array.isArray(responses.options)) {
          options = responses.options;
          searchResult = `${response}: \n`;
        } else {
          searchResult = `${response}: \n ${options}`;
        }
      } else if (target_flie_path.length > 0) {
        // 基因序列文件
        const target_flie_path = responses.target_flie_path;
        const fileList = target_flie_path.map((p) => {
          const url = path.join('/v1/files/xsearchdev', p);
          return `[${path.basename(url)}](${url})`;
        });
        searchResult = `${response}\n${fileList.join('\n')}\n`;
      } else {
        searchResult = response;
      }
      this.logger.debug('searchResult\n' + searchResult);
    } catch (e) {
      this.logger.error(`ncbiSearch error: ${e}`, e.stack);
      searchResult = e.message;
      success = false;
    }

    // yield* this.mockGenerating(`${success ? 'Success' : 'Fail'}\n`);

    yield {
      optionInfo: {
        options,
        cancer_list,
        data_list,
        selected_options: [],
        stage,
        state,
        search_type,
        upload_file,
        species_identification_dict,
      },
      operations,
      content: searchResult,
      role: Role.Assistant,
    };
  }
}
