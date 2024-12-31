import * as WebSocket from 'ws';
import { Logger } from '@nestjs/common';
import { AgentMessageEntity } from '../../entities/agent-message.entity';
import { CacheService } from 'src/cache/cache.service';
import { AgentType } from '@xpcr/common';
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
  data?: any;
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
  data?: any;
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
// submit SpeciesIdentificationInfo
export interface SpeciesIdentificationInfo {
  stage: number;
  search_type?: NCBIFunctionType;
  upload_file_flag?: boolean;
  species_identification_dict?: Record<string, string>;
  operations?: Operation[];
  data?: object;
}
// SpeciesIdentificationResult response
export interface SpeciesIdentificationResult {
  response: string;
  upload_file: boolean;
  stage: number;
  state: 'continue' | 'stop';
  search_type: NCBIFunctionType;
  species_identification_dict: Record<string, string>;
  operations?: Operation[];
  strain_path?: string[];
  non_target_path?: string[];
  target_gene_path?: string[];
  data?: object;
}
// submit PathogenDrugInfo
export interface PathogenDrugInfo {
  stage: number;
  search_type: NCBIFunctionType;
  pathogen_drug_resistance_dict?: Record<string, string>;
  operations?: Operation[];
  data?: object;
}
// PathogenDrugResult response
export interface PathogenDrugResult {
  response: string;
  stage: number;
  state: 'continue' | 'stop';
  search_type: NCBIFunctionType;
  filter_result_path?: string;
  pathogen_drug_resistance_dict?: Record<string, string>;
  operations?: Operation[];
  sequences?: string;
  data?: object;
}

// submit ProteinMutationInfo
export interface ProteinMutationInfo {
  stage: number;
  search_type: NCBIFunctionType;
  protein_mutation_dict?: Record<string, string>;
  operations?: Operation[];
  data?: object;
}
// ProteinMutationResult response
export interface ProteinMutationResult {
  response: string;
  stage: number;
  state: 'continue' | 'stop';
  search_type: NCBIFunctionType;
  protein_mutation_dict?: Record<string, string>;
  protein_gene_seq_path?: string[];
  operations?: Operation[];
  data?: object;
}
// submit CancerOptionInfo
export interface CancerOptionInfo {
  response: string;
  cancer_list: Record<string, any>;
  operations?: Operation[];
  stage: number;
  search_type: NCBIFunctionType;
  search_system_prompt: string;
  state: 'continue' | 'stop';
  selected_options?: Array<string>;
  data?: object;
}
// submit GeneticDisorderInfo
export interface GeneticDisorderInfo {
  response: string;
  operations?: Operation[];
  stage: number;
  search_type: NCBIFunctionType;
  state: 'continue' | 'stop';
  search_system_prompt: string;
  selected_options?: Array<string>;
  data?: object;
}
export interface SnpPrimerDesignInfo {
  operations?: Operation[];
  primer_design_dict: Record<string, string>;
  primer_design_prompt: string;
  primer_type: PrimerDesignFunctionType | string;
  stage: number;
  state: 'continue' | 'stop';
  data?: object;
}

export interface RedesignInfo {
  operations: Operation[];
  primer_type: PrimerDesignFunctionType | string;
  stage: number;
  state: 'continue' | 'stop';
  data?: object;
}

export interface DownLoadTypeOptionsInfo {
  response: string;
  options: Array<Record<string, any>>;
  stage: number;
  data_list: Array<Record<string, any>>;
  search_type: NCBIFunctionType;
  search_system_prompt: string;
  state: 'continue' | 'stop';
  data?: object;
}

export interface Operation {
  title: string;
  key: string;
  options: any;
  value: any[];
  type: string[];
}

// search response
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
  whole_genome_type = 'whole_genome_type',
}
export enum PrimerDesignFunctionType {
  snp_primer_design_type = 'snp_primer_design_type',
  redesign_primer_type = 'redesign_primer_type',
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
  protected domain: string;
  protected resolve?: (value: any) => void;
  protected ws: WebSocket;

  protected connected = false;

  constructor(
    readonly cacheService: CacheService,
    readonly conversationUUID: string,
  ) {
    this.cacheService.setObject(
      `${AgentType.SEQUENCE_SEARCH}:${this.conversationUUID}`,
      false,
    );
    this.initConnect();
  }

  private initConnect() {
    this.ip = process.env.NCBI_HOST;
    this.port = process.env.NCBI_PORT;
    this.path = process.env.NCBI_PATH;
    this.domain = process.env.NCBI_DOMAIN;
    let url = '';
    if (this.domain) {
      url = `wss://${this.domain}${this.path}`;
    } else {
      url = `ws://${this.ip}:${this.port}${this.path}`;
    }
    this.logger.debug('ncbi-search===>initConnect====>', url);
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
      this.cacheService.setObject(
        `${AgentType.SEQUENCE_SEARCH}:${this.conversationUUID}`,
        false,
      );
    });

    this.ws.on('close', () => {
      this.logger.debug('== NCBI Search connect Close ==');
      this.cacheService.setObject(
        `${AgentType.SEQUENCE_SEARCH}:${this.conversationUUID}`,
        false,
      );
    });

    this.ws.on('message', (data) => {
      const text = data.toString();
      if (text === 'heartbeat') {
        this.logger.debug(`== NCBI Search heartbeat ==`);
      } else {
        try {
          this.cacheService.setObject(
            `${AgentType.SEQUENCE_SEARCH}:${this.conversationUUID}`,
            false,
          );
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

  async call(input: string, history: AgentMessageEntity[]) {
    const generating = await this.cacheService.getObject(
      `${AgentType.SEQUENCE_SEARCH}:${this.conversationUUID}`,
    );
    if (generating) {
      return Promise.resolve<NCBIResultContent>({
        responses: {
          response: 'NCBI search is running, please wait...',
        },
      } as NCBIResultContent);
    }
    this.cacheService.setObject(
      `${AgentType.SEQUENCE_SEARCH}:${this.conversationUUID}`,
      true,
    );
    const sendObj = {
      instruction: input,
      history,
      type: 0, // deprecated
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
