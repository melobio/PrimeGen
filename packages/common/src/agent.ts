import {NodeTypes} from './node';
import type { Node } from './node';
export enum AgentType {
  FAULT = 'Fault',
  // Internet Search
  INTERNET_SEARCH = 'Internet Search',
  // Sequence Search
  SEQUENCE_SEARCH = 'Sequence Search',
  PRIMER_DESIGN = 'Primer Design',
  PROTOCOL_DESIGN = 'Protocol Design',
  // Documents Search
  DOCUMENTS_SEARCH = 'Documents Search',
  // Protocol Optimization
  PROTOCOL_OPTIMIZATION = 'Protocol Optimization',
  CODE_EXECUTION = 'Code Execution',
  // Experiment Designer
  EXPERIMENT_DESIGNER = 'Experiment Designer',
  // Active Learning
  ACTIVE_LEARNING = 'Active Learning',
}
// Call By RouterAgent
export const AgentFunctionByType = {
  [AgentType.FAULT]: 'ot2_fault_check',
  [AgentType.CODE_EXECUTION]: 'execute_ot2_protocol',
  [AgentType.INTERNET_SEARCH]: 'internet_search',
  [AgentType.SEQUENCE_SEARCH]: 'sequence_search',
  [AgentType.PRIMER_DESIGN]: 'design_primer',
  [AgentType.PROTOCOL_DESIGN]: 'design_protocol',
}

export const AgentFunctions = {
  CHECK_OT2_STATE: 'check_ot2_state',
  EXECUTE_OT2_PROTOCOL: 'execute_ot2_protocol',
  GOOGLE_SEARCH: 'google_search',
  WIKIPEDIA_SEARCH: 'wikipedia_search',
  NCBI_SEARCH: 'ncbi_search',
  DESIGN_PRIMER: 'design_primer',
  DESIGN_PROTOCOL: 'design_protocol',
}

export interface AgentNode extends Node {
  type: NodeTypes.Agents;
  agentType: AgentType;
}