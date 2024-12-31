import { BaseAgents } from './base-agents';
import { ChatMessage } from '@azure/openai';
import { MiddleMessageType, Role } from '../entities/message.entity';
import { AgentFunctions } from '@xpcr/common';
import {
  NCBIResultContent,
  NcbiSearch,
  TargetGeneInfo,
  NCBIFunctionType,
  CancerOptionInfo,
  GeneticDisorderInfo,
  PathogenDrugInfo,
  SpeciesIdentificationInfo,
  ProteinMutationInfo,
  PathogenDrugResult,
  ProteinMutationResult,
  CancerResult,
} from './ext/ncbi-search';
import { Agents } from '../../tools/agents/entities/agents.entity';
import { DataSource } from 'typeorm';
import { Logger } from '@nestjs/common';
import * as fs from 'fs';
import * as path from 'path';
import { ConversationService } from '../conversation.service';
import { CacheService } from 'src/cache/cache.service';

const createFileListString = (key, responses) => {
  const fileList = responses[key] || [];
  if (fileList.filter((a) => a).length === 0) return '';
  const keyName = key?.replaceAll('_', ' ').toUpperCase() || '';
  const fileListLinks = fileList.map((p) => {
    const url = path.join('/v1/files/xsearchdev', p);
    return `[${path.basename(url)}](${url.replaceAll(' ', '%20')})`;
  });

  return `**${keyName}**: ${fileListLinks.join('\n')};\n`;
};

const createFileString = (key, response) => {
  const file_path = response[key] || '';
  if (file_path) {
    const path_url = path.join('/v1/files/xsearchdev', file_path);
    const path_url_mdStr = `[${path.basename(path_url)}](${path_url.replaceAll(
      ' ',
      '%20',
    )})`;
    return `\n${path_url_mdStr};\n`;
  } else {
    return '';
  }
};

const objectToMarkdownTable = (data) => {
  const keys = Object.keys(data);
  let rows = Object.keys(data[keys[0]]).map((key) => {
    return keys.map((k) => data[k][key] || '-');
  });
  if (rows.length == 0) {
    rows = [keys.map((item) => '-')];
  }
  const header = `| ${keys.join(' | ')} |\n`;
  const divider = `|${keys.map(() => ' :---: ').join('|')}|\n`;
  const body = rows.map((row) => `| ${row.join(' | ')} |`).join('\n');

  return `${header}${divider}${body}`;
};

export class SequenceSearchAgents extends BaseAgents {
  initMessages: ChatMessage[] = [];
  availableFunctions = {
    [AgentFunctions.NCBI_SEARCH]: this.callNCBISearch.bind(this),
  };

  ncbiSearch: NcbiSearch;

  constructor(
    readonly agent: Agents,
    readonly dataSource: DataSource,
    readonly cacheService: CacheService,
    readonly chatsService: ConversationService,
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
    this.ncbiSearch = new NcbiSearch(this.cacheService, this.conversationUUID);
    this.chatsService.setNcbiSearch(this.ncbiSearch, this.conversationUUID);
  }

  // Doesn't need Search Agent Internal LLM to summarize, boost efficiency
  async *send({
    userInput,
    description,
    optionInfo,
  }: {
    userInput: string;
    description?: string;
    optionInfo?:
      | CancerOptionInfo
      | GeneticDisorderInfo
      | PathogenDrugInfo
      | SpeciesIdentificationInfo
      | ProteinMutationInfo;
  }): AsyncGenerator<any, void, any> {
    // yield* this.mockGenerating(`[${this.agent.name}]\n`);
    let content = '';
    let operation = [];
    let newOptionInfo = { state: 'stop' };
    for (const key in this.availableFunctions) {
      const queryFunc = this.availableFunctions[key];
      for await (const msg of queryFunc({
        query: `${userInput}`,
        optionInfo,
      })) {
        if (msg?.optionInfo) {
          newOptionInfo = msg?.optionInfo;
        }
        if (msg?.operation) {
          operation = msg?.operation;
        }
        // Send User Input to Seq Agent
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
      operation,
    };
    await this.saveMessage(`${userInput} ${description}`, Role.User);
    await this.saveMessage(content, Role.Assistant);
  }

  private async *callNCBISearch({
    query,
    optionInfo,
  }: {
    query: string;
    optionInfo?:
      | CancerOptionInfo
      | GeneticDisorderInfo
      | PathogenDrugInfo
      | SpeciesIdentificationInfo
      | ProteinMutationInfo;
  }) {
    this.logger.debug(
      `callNCBISearch==query===> ${query}; optionInfo==> ${JSON.stringify(
        optionInfo,
      )}`,
    );
    let success = true;
    let content = '';
    let responses;
    try {
      const nCBIResult: NCBIResultContent = await this.ncbiSearch.call(
        JSON.stringify(optionInfo),
        this.history,
      );
      if (typeof nCBIResult == 'string') {
        content = `${nCBIResult} \n`;
      } else {
        responses = nCBIResult.responses;
        content = `${responses.response} \n`;
        try {
          if (
            responses?.search_type ==
            NCBIFunctionType.species_identification_type
          ) {
            content += createFileListString(
              'target_cds_path',
              responses.data.species_identification_dict,
            );
          } else if (
            responses?.search_type ==
            NCBIFunctionType.pathogen_drug_resistance_type
          ) {
            const responses: PathogenDrugResult = nCBIResult.responses;
            content += createFileString(`filter_result_path`, responses.data);
            content += createFileString(`sequences`, responses.data);
          } else if (
            responses?.search_type == NCBIFunctionType.protein_mutation_type
          ) {
            const responses: ProteinMutationResult = nCBIResult.responses;
            content += createFileListString(
              'protein_gene_seq_path',
              responses.data,
            );
          } else if (responses?.search_type == NCBIFunctionType.cancer_type) {
            const responses: CancerResult = nCBIResult.responses;
            const target_gene_info = responses?.data?.target_gene_info;
            if (Array.isArray(target_gene_info)) {
              content += createFileListString(
                'target_gene_info',
                responses.data,
              );
            } else if (target_gene_info) {
              const tableStr = objectToMarkdownTable(target_gene_info);
              content += tableStr;
            }
          } else if (
            responses?.search_type == NCBIFunctionType.genetic_disorder_type
          ) {
            const responses: TargetGeneInfo = nCBIResult.responses;
            content += createFileListString('target_gene_path', responses.data);
            const target_gene_info = responses?.data?.target_gene_info;
            if (target_gene_info) {
              const tableStr = objectToMarkdownTable(target_gene_info);
              content += `\n ${tableStr}`;
            }
            if (responses?.data?.history_data?.dna_seqs) {
              content += createFileListString(
                'dna_seqs',
                responses?.data?.history_data,
              );
            }
          } else if (
            responses?.search_type == NCBIFunctionType.whole_genome_type
          ) {
            const responses = nCBIResult.responses;
            content += createFileListString('target_fna_list', responses.data);
          }
          // else if (
          //   responses?.search_type == NCBIFunctionType.download_type
          // ) {
          //   const responses: TargetProteinMutationInfo = nCBIResult.responses;
          //   const protein_gene_seq_path = responses?.protein_gene_seq_path || [];
          //   const fileList = protein_gene_seq_path.map((p) => {
          //     const url = path.join('/v1/files/xsearchdev', p);
          //     return `[${path.basename(url)}](${url})`;
          //   });
          //   agentMessage.content = `${responses.response}\n${fileList.join(
          //     '\n',
          //   )}\n`;
          // }
        } catch (error) {
          this.logger.debug(error);
        }
        if (Array.isArray(responses?.operations)) {
          responses.operations = responses?.operations.filter(
            (item) => typeof item == 'object',
          );
        }
      }
    } catch (e) {
      this.logger.error(`ncbiSearch error: ${e}`, e.stack);
      success = false;
    }
    // yield* this.mockGenerating(`${success ? 'Success' : 'Fail'}\n`);
    yield {
      optionInfo: responses,
      content,
      role: Role.Assistant,
    };
  }
}
