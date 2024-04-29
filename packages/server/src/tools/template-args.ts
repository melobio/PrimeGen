export enum TemplateArgType {
  TextField = 'TextField',
  TextArea = 'TextArea',
  Select = 'Select',
  Secret = 'Secret',
}
export class TemplateArgs {
  label: string;
  key: string;
  value: string;
  type: TemplateArgType;
  options?: string[];
  hint?: string;
  desc?: string;
  constructor(
    label: string,
    key: string,
    value: string,
    type: TemplateArgType,
    options?: string[],
    hint?: string,
    desc?: string,
  ) {
    this.label = label;
    this.key = key;
    this.value = value;
    this.type = type;
    this.options = options;
    this.hint = hint;
    this.desc = desc;
  }
}