import type { RunTimeCommand } from '@xpcr/shared';

export interface CreateProtocolResponse {
  data: {
    id: string;
    protocolType: string;
    robotType: string;
    metadata: {
      author: string;
      protocolName: string;
      apiLevel: string;
    };
  };
}
export interface CreateRunResponse {
  data: {
    id: string;
    status: string;
    current: boolean;
    protocolId: string;
  };
}

type DistributiveOmit<T, K extends keyof T> = T extends any
  ? Omit<T, K>
  : never;
export type RunCommandSummary = DistributiveOmit<RunTimeCommand, 'result'>;

export interface CommandsData {
  data: RunCommandSummary[];
  meta: GetCommandsParams & { totalLength: number };
  links: CommandsLinks;
}

export interface GetCommandsParams {
  cursor: number | null; // the index of the command at the center of the window
  pageLength: number; // the number of items to include
}

export interface CommandsLinks {
  current: {
    // link to the currently executing command
    href: string;
    meta: {
      runId: string;
      commandId: string;
      key: string;
      createdAt: string;
      index: number;
    };
  };
}

export enum RunActionType {
  PLAY = 'play',
  PAUSE = 'pause',
  STOP = 'stop',
}
export type CheckActionType =
  | 'pickUpTip'
  | 'dropTip'
  | 'open_lid'
  | 'close_lid';

export const RUN_STATUS_IDLE = 'idle' as const;
export const RUN_STATUS_RUNNING = 'running' as const;
export const RUN_STATUS_PAUSE_REQUESTED = 'pause-requested' as const;
export const RUN_STATUS_PAUSED = 'paused';
export const RUN_STATUS_STOP_REQUESTED = 'stop-requested' as const;
export const RUN_STATUS_STOPPED = 'stopped' as const;
export const RUN_STATUS_FAILED = 'failed' as const;
export const RUN_STATUS_FINISHING = 'finishing' as const;
export const RUN_STATUS_SUCCEEDED = 'succeeded' as const;
export const RUN_STATUS_BLOCKED_BY_OPEN_DOOR = 'blocked-by-open-door' as const;

export type RunStatus =
  | typeof RUN_STATUS_IDLE
  | typeof RUN_STATUS_RUNNING
  | typeof RUN_STATUS_PAUSE_REQUESTED
  | typeof RUN_STATUS_PAUSED
  | typeof RUN_STATUS_STOP_REQUESTED
  | typeof RUN_STATUS_STOPPED
  | typeof RUN_STATUS_FAILED
  | typeof RUN_STATUS_FINISHING
  | typeof RUN_STATUS_SUCCEEDED
  | typeof RUN_STATUS_BLOCKED_BY_OPEN_DOOR;

export const FINISHED_RUN_STATUSES = [
  RUN_STATUS_SUCCEEDED,
  RUN_STATUS_STOPPED,
  RUN_STATUS_FAILED,
];
