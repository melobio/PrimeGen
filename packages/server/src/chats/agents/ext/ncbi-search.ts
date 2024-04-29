import * as WebSocket from 'ws';
import { Logger } from '@nestjs/common';
import { AgentMessageEntity } from '../../entities/agent-message.entity';
export interface TargetGeneInfo {
  response: string;
  target_gene_info: {
    chrom: { [key: number | string]: number | string };
    start: { [key: number | string]: number | string };
    gene: { [key: number | string]: number | string };
    type: { [key: number | string]: number | string };
    HGNC: { [key: number | string]: number | string };
    end: { [key: number | string]: number | string };
  };
  search_type: NCBIFunctionType;
  search_system_prompt: string;
  state: 'continue' | 'stop';
}

export interface CancerResult {
  response: string;
  target_gene_info: {
    'Gene name': { [key: number | string]: number | string };
    Chromosome: { [key: number | string]: number | string };
    'Gene start': { [key: number | string]: number | string };
    'Gene end': { [key: number | string]: number | string };
    'Start position in gene': { [key: number | string]: number | string };
    'End position in gene': { [key: number | string]: number | string };
  };
  search_type: NCBIFunctionType;
  search_system_prompt: string;
  state: 'continue' | 'stop';
  stage?: number;
}

export interface UploadFileInfo {
  response: string;
  options: Array<Record<string, any>>;
  search_type: NCBIFunctionType;
  search_system_prompt: string;
  state: 'continue' | 'stop';
  upload?: 'True';
  upload_file?: Array<string>;
  upload_title?: string;
}
//物种鉴定提交的数据
export interface SpeciesIdentificationInfo {
  stage: number;
  search_type?: NCBIFunctionType;
  upload_file_flag?: boolean;
  species_identification_dict?: Record<string, string>;
  operations?: Operation[];
}
//物种鉴定的结果
export interface SpeciesIdentificationResult {
  response: string;
  upload_file: boolean;
  stage: number;
  state: 'continue' | 'stop';
  search_type: NCBIFunctionType;
  species_identification_dict: Record<string, string>;
  strain_path?: string;
  operations?: Operation[];
  non_target_path?: string[];
}
// 耐药提交的数据
export interface PathogenDrugInfo {
  stage: number;
  search_type: NCBIFunctionType;
  pathogen_drug_resistance_dict?: Record<string, string>;
  operations?: Operation[];
}
// 耐药的结果
export interface PathogenDrugResult {
  response: string;
  stage: number;
  state: 'continue' | 'stop';
  search_type: NCBIFunctionType;
  filter_result_path?: string;
  pathogen_drug_resistance_dict?: Record<string, string>;
  operations?: Operation[];
  sequences?: string;
}

// 蛋白质突变提交的数据
export interface ProteinMutationInfo {
  stage: number;
  search_type: NCBIFunctionType;
  protein_mutation_dict?: Record<string, string>;
  operations?: Operation[];
}
// 蛋白质突变的结果
export interface ProteinMutationResult {
  response: string;
  stage: number;
  state: 'continue' | 'stop';
  search_type: NCBIFunctionType;
  protein_mutation_dict?: Record<string, string>;
  protein_gene_seq_path?: string[];
  operations?: Operation[];
}
//癌症提交的数据
export interface CancerOptionInfo {
  response: string;
  cancer_list: Record<string, any>;
  stage: number;
  search_type: NCBIFunctionType;
  search_system_prompt: string;
  state: 'continue' | 'stop';
  selected_options?: Array<string>;
}
// 遗传病提交的数据
export interface GeneticDisorderInfo {
  response: string;
  options: Array<Record<string, any>>;
  stage: number;
  search_type: NCBIFunctionType;
  state: 'continue' | 'stop';
  search_system_prompt: string;
  selected_options?: Array<string>;
}

export interface DownLoadTypeOptionsInfo {
  response: string;
  options: Array<Record<string, any>>;
  stage: number;
  data_list: Array<Record<string, any>>;
  search_type: NCBIFunctionType;
  search_system_prompt: string;
  state: 'continue' | 'stop';
}

type ExperimentStrain = 'Group' | 'Microbe';
export interface Operation {
  title: string;
  key: string;
  options: any;
  value: any[];
  type: string[];
}

// search 响应结果
export interface NCBIResultContent {
  responses: ProteinMutationInfo &
    DownLoadTypeOptionsInfo &
    GeneticDisorderInfo &
    CancerOptionInfo &
    PathogenDrugResult &
    TargetGeneInfo &
    CancerResult &
    UploadFileInfo & {
      target_flie_path: string[];
    } & SpeciesIdentificationResult &
    ProteinMutationResult;
  history_conversation: Array<any>[];
}
export enum NCBIFunctionType {
  cancer_type = 'cancer_type',
  genetic_disorder_type = 'genetic_disorder_type',
  protein_mutation_type = 'protein_mutation_type',
  species_identification_type = 'species_identification_type',
  pathogen_drug_resistance_type = 'pathogen_drug_resistance_type',
  download_type = 'download_type',
}

export type NCBIStateType = {
  continue: 'continue';
  stop: 'stop';
};
export interface NCBIResult {
  content: NCBIResultContent;
  type: number;
  state?: string;
}

export class NcbiSearch {
  logger = new Logger(NcbiSearch.name);
  protected ip: string;
  protected port: string;
  protected path: string;

  protected generating = false;
  protected generatingContent = '';
  protected resolve?: (value: any) => void;
  protected ws: WebSocket;

  protected connected = false;

  constructor() {
    this.initConnect();
  }

  private initConnect() {
    this.ip = process.env.NCBI_HOST;
    this.port = process.env.NCBI_PORT;
    this.path = process.env.NCBI_PATH;

    const url = `ws://${this.ip}:${this.port}${this.path}`;
    this.ws = new WebSocket(url);

    this.ws.on('open', () => {
      this.logger.log(`== NCBI Search connect success. ==`);
      this.connected = true;
    });

    this.ws.on('error', (error) => {
      this.logger.error('== NCBI Search connect error ==');
      this.logger.error(`${error.toString()}`);
      this.sendResolve({
        responses: {
          response: `NCBI Search error: ${error.toString()}`,
        },
      } as NCBIResultContent);
      this.generating = false;
    });

    this.ws.on('close', () => {
      this.logger.debug('== NCBI Search connect Close ==');
      this.generating = false;
    });

    this.ws.on('message', (data) => {
      const text = data.toString();
      if (text === 'heartbeat') {
        this.logger.debug(`== NCBI Search heartbeat ==`);
      } else {
        try {
          this.generating = false;
          const ncbiResult: NCBIResult = JSON.parse(text);
          this.logger.debug('== NCBI Search Message ==');
          this.logger.debug(JSON.parse(text));
          if (typeof ncbiResult.content == 'string') {
            const errorTip = `We apologize to the failure to help immediately;It seems that our assistants have encountered a little problem;Please try to check if there is an unsubmitted operation, our assistant will be happy to serve you.`;
            this.sendResolve(errorTip);
          } else {
            const state: 'stop' | 'continue' =
              ncbiResult.content?.responses?.state;
            if (state === 'stop' || state === 'continue') {
              this.sendResolve(ncbiResult.content);
            }
          }
        } catch (e) {
          this.logger.debug(`== NCBI Search handle message error ==`);
          this.logger.debug(e);
        }
      }
    });
    this.logger.log(`NCBI search Init Success`);
  }
  private sendResolve(content: NCBIResultContent | string) {
    if (this.resolve) {
      this.resolve(content);
      this.resolve = null;
    }
  }
  call(input: string, history: AgentMessageEntity[]) {
    if (this.generating) {
      return Promise.resolve<NCBIResultContent>({
        responses: {
          response: 'NCBI search is running, please wait...',
        },
      } as NCBIResultContent);
    }
    this.generating = true;
    this.generatingContent = '';
    const sendObj = {
      instruction: input,
      history,
      type: 0, // 无用字段
    };
    this.ws.send(JSON.stringify(sendObj));
    this.logger.log(`send: ${JSON.stringify(sendObj)}`);

    return new Promise<NCBIResultContent>((resolve) => {
      this.resolve = resolve;
    });
  }

  callFromService(input: string, history: AgentMessageEntity[]) {
    if (this.generating) {
      return Promise.resolve<NCBIResultContent>({
        responses: {
          response: 'NCBI search is running, please wait...',
        },
      } as NCBIResultContent);
    }
    this.generating = true;
    this.generatingContent = '';
    const sendObj = {
      instruction: input,
      history,
      type: 0, // 无用字段
    };
    if (!this.ws) {
      this.initConnect();
    }
    this.ws.send(JSON.stringify(sendObj));
    this.logger.log(`send: ${JSON.stringify(sendObj)}`);

    return new Promise<NCBIResultContent>((resolve) => {
      this.resolve = resolve;
    });
  }
}
