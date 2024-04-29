import type { NodeChild } from "./node";

export enum ExperimentsState {
  CREATED = 'CREATED',
  RUNNING = 'RUNNING',
  PAUSED = 'PAUSED',
  FINISHED = 'FINISHED',
  CANCELED = 'CANCELED',
  FAILED = 'FAILED',
}

export interface ExperimentInterface {
  id: number;
  uuid: string;
  name: string;
  nodeChildren: NodeChild[];
  desc: string;
  creator: string;
  creatorName?: string;
  state: ExperimentsState;
  createTime: Date;
  updateTime: Date;
}