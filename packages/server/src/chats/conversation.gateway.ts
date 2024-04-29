import {
  OnGatewayConnection,
  OnGatewayDisconnect,
  SubscribeMessage,
  WebSocketGateway,
  WebSocketServer,
} from '@nestjs/websockets';
import { Server, Socket } from 'socket.io';
import { Logger, OnModuleInit } from '@nestjs/common';
import { ConversationService } from './conversation.service';
import { Role } from './entities/message.entity';
import {
  JETSON_BAD_TIP,
  TEXT_ANSWER_DONE,
  TEXT_QUESTION,
  VOICE_QUESTION,
  SEND_OPTION_NCBI_SEARCH,
} from './conversation.constant';
import axios, { AxiosRequestConfig } from 'axios';
import type { AxiosInstance } from 'axios';
import * as FormData from 'form-data';
import * as fs from 'fs';
import {
  CancerOptionInfo,
  GeneticDisorderInfo,
  PathogenDrugInfo,
  SpeciesIdentificationInfo,
  ProteinMutationInfo,
} from './agents/ext/ncbi-search';

@WebSocketGateway({
  path: '/pcr-ws',
  allowEIO3: true,
  transports: ['polling', 'websocket'],
  cors: true,
})
export class ConversationGateway
  implements OnGatewayConnection, OnGatewayDisconnect, OnModuleInit
{
  logger = new Logger(ConversationGateway.name);
  @WebSocketServer()
  server: Server;

  whisperClient: AxiosInstance;

  constructor(private readonly chatsService: ConversationService) {
    this.whisperClient = axios.create({
      timeout: 40000,
    });
  }

  handleDisconnect(client: Socket) {
    this.logger.log('handleDisconnect ' + client.id);
  }

  handleConnection(client: Socket) {
    this.chatsService.setClient(client);
    this.logger.log('handleConnection ' + client.id);
  }

  onModuleInit() {
    this.logger.log('onModuleInit');
  }

  private async voiceToText(voice: Buffer) {
    const data = new FormData();
    data.append('file', voice, { filename: 'test.webm' });
    data.append('language', 'zh');
    const url = `${process.env.WHISPER_AZURE_OPENAI_ENDPOINT}/openai/deployments/${process.env.WHISPER_DEPLOYMENT_NAME}/audio/transcriptions?api-version=2023-09-01-preview`;
    const config: AxiosRequestConfig = {
      method: 'post',
      maxBodyLength: Infinity,
      url,
      headers: {
        'api-key': process.env.WHISPER_AZURE_OPENAI_KEY,
        ...data.getHeaders(),
      },
      data: data,
    };
    try {
      const response = await axios.request(config);
      this.logger.log(response.data);
      return response.data.text;
    } catch (e) {
      console.error(e);
      return '';
    }
  }
  @SubscribeMessage(VOICE_QUESTION)
  async voiceQuestion(
    client: Socket,
    { voice, conversationUUID }: { voice: string; conversationUUID: string },
  ) {
    // this.logger.log('voiceQuestion: base64 len: ' + voice);
    const voiceBuffer = Buffer.from(voice, 'base64');
    // fs.writeFileSync('./test.webm', voiceBuffer);
    this.logger.log(`voiceQuestion: ${voiceBuffer.length}`);
    const text = await this.voiceToText(voiceBuffer);

    if (text) {
      return await this.textQuestion(client, { text, conversationUUID });
    }
  }

  @SubscribeMessage(TEXT_QUESTION)
  async textQuestion(
    client: Socket,
    {
      text,
      conversationUUID,
    }: {
      text: string;
      conversationUUID: string;
    },
  ) {
    this.logger.log('inputText: ' + text);
    setTimeout(() => {
      // 发送text到openai
      this.chatsService
        .sendText(conversationUUID, client, text)
        .catch((e) => {
          this.logger.error(`sendText error: ${e.message}`, e.stack);
        })
        .finally(() => {
          this.logger.log(`sendText done.`);
          client.emit(TEXT_ANSWER_DONE, { conversationUUID });
        });
    }, 200);
    // 创建一个user消息
    const userMessage = await this.chatsService.createMessage(
      '',
      text,
      conversationUUID,
      Role.User,
    );
    return { success: true, data: [userMessage] };
  }

  @SubscribeMessage(JETSON_BAD_TIP)
  async handleBadTip(
    client: Socket,
    { runId, option }: { runId: string; option: 'play' | 'stop' },
  ) {
    const res = await this.chatsService.handleBadTip(option, runId);
    return res;
  }

  @SubscribeMessage(SEND_OPTION_NCBI_SEARCH)
  async handleGeneticDiseasesOptionsSubmit(
    client: Socket,
    {
      text,
      conversationUUID,
      option,
      messageId,
    }: {
      text: string;
      conversationUUID: string;
      option?:
        | CancerOptionInfo
        | GeneticDisorderInfo
        | PathogenDrugInfo
        | SpeciesIdentificationInfo
        | ProteinMutationInfo;
      messageId: number;
    },
  ) {
    setTimeout(async () => {
      await this.chatsService.updateNcbiSearchMessage({
        messageId,
        option,
      });
      const res = await this.chatsService.handleNcbiSearch({
        option,
        conversationUUID,
      });
      await this.chatsService.updateNcbiSearchMessage({
        messageId,
        option,
        optionRes: res.data,
      });
      if (res.data.state == 'stop' && res.success) {
        const routerAgents = await this.chatsService.getOrCreateRouterAgents(
          conversationUUID,
        );
        await routerAgents.handleNcbiSearchStop(client, res.data);
      }
      client.emit(TEXT_ANSWER_DONE, {
        conversationUUID,
      });
    }, 200);
    const userMessage = await this.chatsService.createMessage(
      '',
      text,
      conversationUUID,
      Role.User,
    );
    return { success: true, data: [userMessage] };
  }
}
