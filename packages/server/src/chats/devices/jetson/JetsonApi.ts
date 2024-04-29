import axios, { AxiosRequestConfig } from 'axios';
import type { AxiosInstance } from 'axios';
import { Logger } from '@nestjs/common';
import * as WebSocket from 'ws';

export interface JetsonApiResponse<T> {
  data: T;
  message: string;
  success: boolean;
}

export class JetsonApi {
  ws: WebSocket;
  logger = new Logger(JetsonApi.name);
  client: AxiosInstance;
  baseURL: string;
  constructor(baseURL: string) {
    this.baseURL = baseURL;
    this.client = axios.create({
      baseURL,
      timeout: 40000,
    });
    this.client.interceptors.request.use((request) => {
      return request;
    });
    this.logger.log(`init jetson client success: ${baseURL}`);
  }

  wsConnectInit() {
    const ip = process.env.JETSON_API_BASE_HOST;
    const port = process.env.JETSON_API_BASE_PORT;
    const path = process.env.JETSON_API_BASE_PATH;
    const url = `ws://${ip}:${port}${path}`;
    this.ws = new WebSocket(url);
    this.ws.on('open', () => {
      this.logger.log(`connect success.`);
    });
    this.ws.on('error', (error) => {
      this.logger.error(`error: ${error.toString()}`);
    });

    this.ws.on('close', () => {
      this.ws = null;
      this.logger.log(`close`);
    });
  }

  async request<T = any>(config: AxiosRequestConfig) {
    try {
      const { data } = await this.client.request<T>(config);
      return {
        success: true,
        data,
        message: 'success',
      } as JetsonApiResponse<T>;
    } catch (e) {
      return {
        success: false,
        data: null,
        message: e.toString(),
      } as JetsonApiResponse<null>;
    }
  }
  async check<T>(
    checkType: 'tips' | 'PCR' | 'start' | 'stop',
    action: 'pickUpTip' | 'dropTip' | 'open_lid' | 'close_lid',
    single: boolean,
  ) {
    return await this.request<T>({
      url: '/checkType',
      data: {
        checkType,
        action,
        Single: single,
      },
      method: 'POST',
    });
  }
}
