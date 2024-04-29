export enum NodeFieldType {
  TextArea = 'TextArea',
  Inputs = 'Inputs',
  TextField = 'TextField',
  Select = 'Select',
  Secret = 'Secret',
}

export interface NodeField {
  label: string;
  key: string;
  value?: string;
  type: NodeFieldType;
  options?: string[];
  hint?: string;
  desc?: string;
  required: boolean;
  accept: {};
}
