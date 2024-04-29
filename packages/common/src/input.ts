import {Node, NodeTypes} from "./node";

export enum InputType {
  InputNode = 'Input Node',
  InputFileNode = 'Input File Node',
}
export enum InputFields {
  INPUT = 'input',
}
export interface InputNode extends Node {
  type: NodeTypes.Input;
  inputType: InputType;
}