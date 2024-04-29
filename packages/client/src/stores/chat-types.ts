import type { CompletedProtocolAnalysis } from "@xpcr/shared";
export enum NCBIFunctionType  {
  cancer_type = 'cancer_type',
  genetic_disorder_type ='genetic_disorder_type',
  protein_mutation_type= 'protein_mutation_type',
  species_identification_type='species_identification_type',
  pathogen_drug_resistance_type= 'pathogen_drug_resistance_type',
  download_type='download_type'
}
export interface Operation {
  title: string;
  key: string;
  options: any;
  value: any[];
  type: string[];
}
//物种鉴定

export interface SpeciesIdentificationInfo{
  stage: number;
  search_type?: NCBIFunctionType;
  upload_file_flag?:boolean;
  species_identification_dict?: Record<string, string>;
  operations?: Operation[];
}
export interface SpeciesIdentificationResult {
  response: string;
  upload_file: boolean;
  stage: number;
  state: 'continue' | 'stop';
  search_type: NCBIFunctionType;
  species_identification_dict: Record<string, string>;
  strain_path?: string;
  operations?: Operation[];
}
//耐药
export interface PathogenDrugInfo{
  stage: number;
  search_type: NCBIFunctionType;
  pathogen_drug_resistance_dict?:Record<string, string>;
  operations?: Operation[];
}
export interface PathogenDrugResult {
  response: string;
  stage: number;
  state: 'continue' | 'stop';
  search_type: NCBIFunctionType;
  pathogen_drug_resistance_dict?: Record<string, string>;
  operations?: Operation[];
}
//蛋白质突变
export interface ProteinMutationInfo {
  stage: number;
  search_type: NCBIFunctionType;
  protein_mutation_dict?:Record<string, string>;
  operations?: Operation[];
}
//蛋白质突变的结果
export interface ProteinMutationResult {
  response: string;
  stage: number;
  state: 'continue' | 'stop';
  search_type: NCBIFunctionType;
  protein_mutation_dict?: Record<string, string>;
  protein_gene_seq_path?: string[];
  operations?: Operation[];
}

export interface GeneticDisorderInfo{
  stage: number;
  state: 'continue' | 'stop';
  search_type: NCBIFunctionType;
  options: Array<Record<string, any>>;
  selected_options: Array<string>;
  selectedOptionsInfo?:Array<Record<string, string>>;
}
export interface CancerOptionInfo{
  stage: number;
  search_type: NCBIFunctionType;
  cancer_dict?: Record<string, string>;
  operations?: Operation[];
}
export interface CancerResult{
  response: string;
  stage: number;
  state: 'continue' | 'stop';
  search_type: NCBIFunctionType;
  cancer_dict?: Record<string, string>;
  operations?: Operation[];
}
export interface Chat {
  id: number;
  uuid: string;
  name: string;
  desc: string;
  createTime: string;
  updateTime: string;
  isGenerating?:boolean;
  isRecordingAudio?:boolean
}

export interface Message {
  id: number;
  uuid: string;
  conversationUUID: string;
  role: string;
  header: string;
  agentType: string;
  protocolAnalysis?: CompletedProtocolAnalysis;
  currentCommandIndex?: number;
  faultCheckResult?: {
    success: boolean;
    data: string;
    plot: number[];
    commandIndex?: number;
  }[];
  content: string;
  createTime: string;
  updateTime: string;
  optionInfo?:GeneticDisorderInfo & 
  CancerOptionInfo & 
  CancerResult &
  PathogenDrugResult & 
  ProteinMutationResult &
  SpeciesIdentificationResult;
  generating: boolean;
}

export interface ConversationUUIDtoPrimerFileID {
  conversationUUID: string,
  experimentId: string,
}

export const VOICE_QUESTION = 'VOICE_QUESTION';
export const TEXT_QUESTION = 'TEXT_QUESTION';
export const TEXT_ANSWER_GENERATING = 'TEXT_ANSWER_GENERATING';
export const TEXT_ANSWER_START = 'TEXT_ANSWER_START';
export const TEXT_ANSWER_DONE = 'TEXT_ANSWER_DONE';
export const TEXT_ANSWER_ERROR = 'TEXT_ANSWER_ERROR';
export const TEXT_ANSWER_CREATE = 'TEXT_ANSWER_CREATE';
export const JETSON_BAD_TIP = 'JETSON_BAD_TIP'; //出现坏针头，提示用户
export const INITIATIVE_START_PRIMER_DESIGN = 'INITIATIVE_START_PRIMER_DESIGN'; // 提示是否开始引物设计
export const UPLOAD_PRIMER_EXCEL = 'UPLOAD_PRIMER_EXCEL';
export const SEND_OPTION_NCBI_SEARCH = 'SEND_OPTION_NCBI_SEARCH'