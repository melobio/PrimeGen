export interface FunctionArgItem {
  description: string;
  type: string;
}

export interface FunctionItem {
  name: string;
  description: string;
  parameters: {
    type: 'object';
    properties: {
      [name: string]: FunctionArgItem;
    };
  };
  required: string[];
}