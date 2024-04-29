import {Node, NodeTypes} from './node';

export enum ToolType {
  GOOGLE_SEARCH = 'Google Search',
  BING_SEARCH = 'Bing Search',
  WIKIPEDIA_SEARCH = 'Wikipedia Search',
}

export enum ToolFields {
  GOOGLE_API_KEY = 'googleApiKey',
  GOOGLE_ENGINE_ID = 'googleEngineId',
  WIKIPEDIA_TOP_K_RESULTS = 'wikipediaTopKResults',
  WIKIPEDIA_MAX_DOC_CONTENT_LENGTH = 'wikipediaMaxDocContentLength',
}

export interface ToolNode extends Node {
  type: NodeTypes.Tools;
  toolType: ToolType;
}