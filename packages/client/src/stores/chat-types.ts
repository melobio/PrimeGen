import type { CompletedProtocolAnalysis } from "@xpcr/shared";
export enum NCBIFunctionType  {
  cancer_type = 'cancer_type',
  genetic_disorder_type ='genetic_disorder_type',
  protein_mutation_type= 'protein_mutation_type',
  species_identification_type='species_identification_type',
  pathogen_drug_resistance_type= 'pathogen_drug_resistance_type',
  download_type='download_type'
}

export enum ProtocolDesignFunctionType {
  auto_protocol_design = 'auto_protocol_design',
  template_protocol_design = 'template_protocol_design',
}

export enum PrimerDesignFunctionType  {
  snp_primer_design_type = 'snp_primer_design_type',
}

export interface Operation {
  title: string;
  key: string;
  options: any;
  value: any[];
  type: string[];
  widget?: string;
}
//物种鉴定

export interface SpeciesIdentificationInfo{
  stage: number;
  search_type?: NCBIFunctionType;
  upload_file_flag?:boolean;
  species_identification_dict?: Record<string, string>;
  operations?: Operation[];
  data?: object;
}
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
//耐药
export interface PathogenDrugInfo{
  stage: number;
  search_type: NCBIFunctionType;
  pathogen_drug_resistance_dict?:Record<string, string>;
  operations?: Operation[];
  data?: object;
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
  data?: object;
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
  data?: any;
}
//蛋白质突变范围结果
export interface ProteinMutationLengthResult {
  response: string;
  stage: number,
  state: 'continue' | 'stop';
  search_type: NCBIFunctionType;
  max_length?: number;
  protein_gene_seq_path?: string[];
  operations?: Operation[];
  protein_mutation_dict?: Record<string, string>;
  data?: object;
}
// 蛋白质范围选择
export interface ProteinMutationLength {
  stage: number;
  search_type: NCBIFunctionType;
  start_pos: number;
  end_pos: number;
  protein_gene_seq_path?: string[];
  protein_mutation_dict?: Record<string, string>;
  data?: object;
}

// Protocol Design 入参
export interface ProtocolDesignInfo {
  stage: number;
  protocol_design_type: ProtocolDesignFunctionType;
  operations?: Operation[];
  data?: object;
}

// Protocol Design 出参
export interface ProtocolDesignResultInfo {
  response: string;
  stage: number;
  protocol_design_type: ProtocolDesignFunctionType;
  state: 'continue' | 'stop';
  operations?: Operation[];
  data?: object;
}

export interface GeneticDisorderInfo{
  stage: number;
  state: 'continue' | 'stop';
  search_type: NCBIFunctionType;
  // options: Array<Record<string, any>>;
  // selected_options: Array<string>;
  operations?: Operation[];
  selectedOptionsInfo?:Array<Record<string, string>>;
  data?: object;
}

export interface GeneticLengthInfo{
  stage: number;
  search_type: NCBIFunctionType;
  start_pos: number;
  end_pos: number;
  target_gene_path: Array<string>;
  data?: object;
}

export interface CancerOptionInfo{
  stage: number;
  search_type: NCBIFunctionType;
  cancer_dict?: Record<string, string>;
  operations?: Operation[];
  data?: object;
}
export interface CancerResult{
  response: string;
  stage: number;
  state: 'continue' | 'stop';
  search_type: NCBIFunctionType;
  cancer_dict?: Record<string, string>;
  operations?: Operation[];
  data?: object;
}
export interface SnpPrimerDesignInfo{
  operations?: Operation[];
  primer_design_dict: Record<string,string>,
  primer_design_prompt:string;
  primer_type: PrimerDesignFunctionType;
  stage: number;
  state: 'continue' | 'stop';
  submitted?:boolean;
  data?: object;
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
  ProteinMutationLengthResult &
  ProteinMutationLength &
  ProtocolDesignInfo &
  ProtocolDesignResultInfo &
  SpeciesIdentificationResult &
  SnpPrimerDesignInfo;
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
export const UPDATE_CONVERSATIONS = 'UPDATE_CONVERSATIONS';
export const UPDATE_MESSAGE = 'UPDATE_MESSAGE';
export const STOP_GENERATING_MESSAGE = 'STOP_GENERATING_MESSAGE';
export const SUMMARY_MESSAGE = 'SUMMARY_MESSAGE';
export const TEXT_ANSWER_ERROR = 'TEXT_ANSWER_ERROR';
export const TEXT_ANSWER_CREATE = 'TEXT_ANSWER_CREATE';
export const JETSON_BAD_TIP = 'JETSON_BAD_TIP'; //出现坏针头，提示用户
export const INITIATIVE_START_PRIMER_DESIGN = 'INITIATIVE_START_PRIMER_DESIGN'; // 提示是否开始引物设计
export const RESTART_CONVERSATION = 'RESTART_CONVERSATION'; //重新开始会话
export const UPLOAD_PRIMER_EXCEL = 'UPLOAD_PRIMER_EXCEL';
export const SEND_OPTION_NCBI_SEARCH = 'SEND_OPTION_NCBI_SEARCH'
export const SEND_OPTION_PRIMER_DESIGN = 'SEND_OPTION_PRIMER_DESIGN'
export const SEND_OPTION_PROTOCOL_DESIGN = 'SEND_OPTION_PROTOCOL_DESIGN';