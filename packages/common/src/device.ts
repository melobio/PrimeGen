import {Node, NodeTypes} from './node';

export enum DeviceType {
  OT2 = 'Opentrons OT-2',
  MGI_ALPHA_TOOL = 'MGI AlphaTool',
}
export enum DeviceFields {
  API_BASE = 'apiBase',
  //jetson_api_base
  JETSON_API_BASE = 'jetsonApiBase',
}
export interface DeviceNode extends Node {
  type: NodeTypes.Devices;
  deviceType: DeviceType;
}