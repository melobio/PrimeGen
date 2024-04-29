import axios, { AxiosRequestConfig } from 'axios';
import type { AxiosInstance } from 'axios';
import { Logger } from '@nestjs/common';
import * as uuid from 'uuid';
import { CommandsData, RunActionType } from './ot-types';
import * as fs from 'fs';
import * as path from 'path';

export interface OTApiResponse<T> {
  data: T;
  message: string;
  success: boolean;
}

export class OTApi {
  logger = new Logger(OTApi.name);
  client: AxiosInstance;
  constructor(private readonly baseURL: string) {
    this.client = axios.create({
      baseURL,
      timeout: 40000,
    });
    this.client.interceptors.request.use((request) => {
      request.headers['Opentrons-Version'] = '3';
      return request;
    });
    this.logger.log(`init opentrons client success: ${baseURL}`);
  }
  async request<T = any>(config: AxiosRequestConfig) {
    try {
      const { data } = await this.client.request<T>(config);
      return { success: true, data, message: 'success' } as OTApiResponse<T>;
    } catch (e) {
      return {
        success: false,
        data: null,
        message: JSON.stringify(e.response?.data ?? e.toString()),
      } as OTApiResponse<null>;
    }
  }

  health() {
    return this.client.get<any>(`/health`);
  }

  createProtocol(protocol: string, protocolKey?: string) {
    const formData = new FormData();
    console.log('protocol', protocol);
    formData.append('files', new Blob([protocol]), `${uuid.v4()}.py`);
    if (protocolKey != null) formData.append('key', protocolKey);

    let pendingResolve: (res: { success: boolean; data: string }) => void;
    let pendingReject: (reason: any) => void;
    const promise = new Promise<{ success: boolean; data: string }>(
      (resolve, reject) => {
        pendingResolve = resolve;
        pendingReject = reject;
      },
    );

    this.client
      .post(`/protocols`, formData)
      .then((r) => {
        // console.log(r.data);
        pendingResolve({ success: true, data: JSON.stringify(r.data) });
      })
      .catch((e) => {
        pendingResolve({
          success: false,
          data: JSON.stringify(e.response?.data ?? e.toString()),
        });
      });

    return promise;
  }

  createRun(protocolId: string) {
    const labwareOffsets = fs.readFileSync(
      path.join(process.cwd(), 'assets', 'ot2', 'labwareOffsets.json'),
    );

    const createRunData = {
      labwareOffsets: JSON.parse(String(labwareOffsets)),
      protocolId,
    };
    let pendingResolve: (res: { success: boolean; data: string }) => void;
    let pendingReject: (reason: any) => void;
    const promise = new Promise<{ success: boolean; data: string }>(
      (resolve, reject) => {
        pendingResolve = resolve;
        pendingReject = reject;
      },
    );
    this.client
      .post(`/runs`, { data: createRunData })
      .then((r) => {
        // console.log(r.data);
        pendingResolve({ success: true, data: JSON.stringify(r.data) });
      })
      .catch((e) => {
        pendingResolve({
          success: false,
          data: JSON.stringify(e.response?.data ?? e.toString()),
        });
      });

    return promise;
  }

  async getRun(runId: string) {
    return await this.request({
      url: `/runs/${runId}`,
    });
  }

  runAction(runId: string, action: RunActionType) {
    const createRunActionData = {
      actionType: action,
    };
    let pendingResolve: (res: { success: boolean; data: string }) => void;
    let pendingReject: (reason: any) => void;
    const promise = new Promise<{ success: boolean; data: string }>(
      (resolve, reject) => {
        pendingResolve = resolve;
        pendingReject = reject;
      },
    );

    this.client
      .post(`/runs/${runId}/actions`, { data: createRunActionData })
      .then((r) => {
        // console.log(r.data);
        pendingResolve({ success: true, data: JSON.stringify(r.data) });
      })
      .catch((e) => {
        pendingResolve({
          success: false,
          data: JSON.stringify(e.response?.data ?? e.toString()),
        });
      });

    return promise;
  }

  getRunCommands(runId: string, pageLength: number) {
    // /runs/{runId}/commands
    let pendingResolve: (res: { success: boolean; data: CommandsData }) => void;
    let pendingReject: (reason: any) => void;
    const promise = new Promise<{ success: boolean; data: CommandsData }>(
      (resolve, reject) => {
        pendingResolve = resolve;
        pendingReject = reject;
      },
    );

    this.client
      .get<CommandsData>(`/runs/${runId}/commands`, { params: { pageLength } })
      .then((r) => {
        // console.log(r.data);
        pendingResolve({ success: true, data: r.data });
      })
      .catch((e) => {
        pendingResolve({
          success: false,
          data: e.response?.data ?? e.toString(),
        });
      });

    return promise;
  }

  getProtocolAnalysis(protocolId: string) {
    let pendingResolve: (res: { success: boolean; data: string }) => void;
    let pendingReject: (reason: any) => void;
    const promise = new Promise<{ success: boolean; data: string }>(
      (resolve, reject) => {
        pendingResolve = resolve;
        pendingReject = reject;
      },
    );

    this.client
      .get(`/protocols/${protocolId}/analyses`, {})
      .then((r) => {
        // console.log(r.data);
        pendingResolve({ success: true, data: JSON.stringify(r.data) });
      })
      .catch((e) => {
        pendingResolve({
          success: false,
          data: JSON.stringify(e.response?.data ?? e.toString()),
        });
      });

    return promise;
  }
}
