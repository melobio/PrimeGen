import {
  OnGatewayConnection,
  OnGatewayDisconnect,
  SubscribeMessage,
  WebSocketGateway,
  WebSocketServer,
  WsException,
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
  RESTART_CONVERSATION,
  STOP_GENERATING_MESSAGE,
  TEXT_ANSWER_GENERATING,
  SEND_OPTION_PRIMER_DESIGN,
  SEND_OPTION_PROTOCOL_DESIGN,
} from './conversation.constant';
import axios, { AxiosRequestConfig } from 'axios';
import type { AxiosInstance } from 'axios';
import * as FormData from 'form-data';
import {
  CancerOptionInfo,
  GeneticDisorderInfo,
  PathogenDrugInfo,
  SpeciesIdentificationInfo,
  ProteinMutationInfo,
  SnpPrimerDesignInfo,
} from './agents/ext/ncbi-search';
import { ProtocolDesignInfo } from './agents/ext/protocol-design';
import { AgentType } from '@xpcr/common';
import * as pako from 'pako';

@WebSocketGateway({
  path: '/xpcr-ws',
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
    // 创建一个user消息
    const userMessage = await this.chatsService.createMessage(
      '',
      text,
      conversationUUID,
      Role.User,
    );
    this.chatsService.setCacheObject(
      `${TEXT_ANSWER_GENERATING}:${conversationUUID}`,
      true,
    );
    // 发送text到openai
    this.chatsService
      .sendText(conversationUUID, client, text)
      .catch((e) => {
        this.logger.error(`sendText error: ${e.message}`, e.stack);
      })
      .finally(() => {
        this.chatsService.emitToCliet(client, TEXT_ANSWER_DONE, {
          conversationUUID,
        });
      });
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

  @SubscribeMessage(STOP_GENERATING_MESSAGE)
  async handleStopGeneratingMessage(
    client: Socket,
    { conversationUUID }: { conversationUUID: string },
  ) {
    this.logger.debug('STOP_GENERATING_MESSAGE');
    this.chatsService.resetCacheStatus(conversationUUID);
    this.chatsService.setCacheObject(
      `${TEXT_ANSWER_GENERATING}:${conversationUUID}`,
      false,
    );
    // todo: 停止生成信息后，拦截数据库更新数据
    return { data: [], success: true };
  }

  @SubscribeMessage(SEND_OPTION_NCBI_SEARCH)
  async handleNCBISearchOption(client: Socket, compressedData: Uint8Array) {
    const {
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
    } = JSON.parse(pako.inflate(compressedData, { to: 'string' }));
    this.chatsService.setCacheObject(
      `${TEXT_ANSWER_GENERATING}:${conversationUUID}`,
      true,
    );
    const userMessage = await this.chatsService.createMessage(
      '',
      text,
      conversationUUID,
      Role.User,
    );
    setTimeout(async () => {
      await this.chatsService.updateAgentMessage({
        messageId,
        agentType: AgentType.SEQUENCE_SEARCH,
        option,
      });
      const routerAgents = await this.chatsService.getOrCreateRouterAgents(
        conversationUUID,
      );
      await routerAgents.sendToSearch(client, text, option);
      client.emit(TEXT_ANSWER_DONE, {
        conversationUUID,
      });
    }, 200);
    return { success: true, data: [userMessage] };
  }

  @SubscribeMessage(SEND_OPTION_PRIMER_DESIGN)
  async handlePrimerDesignOption(
    client: Socket,
    {
      text,
      conversationUUID,
      option,
      messageId,
    }: {
      text: string;
      conversationUUID: string;
      option?: SnpPrimerDesignInfo;
      messageId: number;
    },
  ) {
    this.chatsService.setCacheObject(
      `${TEXT_ANSWER_GENERATING}:${conversationUUID}`,
      true,
    );
    const userMessage = await this.chatsService.createMessage(
      '',
      text,
      conversationUUID,
      Role.User,
    );
    setTimeout(async () => {
      await this.chatsService.updateAgentMessage({
        messageId,
        agentType: AgentType.PRIMER_DESIGN,
        option,
      });
      const routerAgents = await this.chatsService.getOrCreateRouterAgents(
        conversationUUID,
      );
      await routerAgents.sendToPrimerDesign(client, text, option);
      client.emit(TEXT_ANSWER_DONE, {
        conversationUUID,
      });
    }, 200);
    return { success: true, data: [userMessage] };
  }

  @SubscribeMessage(SEND_OPTION_PROTOCOL_DESIGN)
  async handleProtocolDesignOption(
    client: Socket,
    {
      text,
      conversationUUID,
      option,
      messageId,
    }: {
      text: string;
      conversationUUID: string;
      option?: ProtocolDesignInfo;
      messageId: number;
    },
  ) {
    this.chatsService.setCacheObject(
      `${TEXT_ANSWER_GENERATING}:${conversationUUID}`,
      true,
    );
    const userMessage = await this.chatsService.createMessage(
      '',
      text,
      conversationUUID,
      Role.User,
    );
    setTimeout(async () => {
      await this.chatsService.updateAgentMessage({
        messageId,
        agentType: AgentType.PROTOCOL_DESIGN,
        option,
      });
      const routerAgents = await this.chatsService.getOrCreateRouterAgents(
        conversationUUID,
      );
      await routerAgents.sendToProtocolDesign(client, text, option);
      client.emit(TEXT_ANSWER_DONE, {
        conversationUUID,
      });
    }, 200);
    return { success: true, data: [userMessage] };
  }

  @SubscribeMessage(RESTART_CONVERSATION)
  async handleRestartConversation(
    client: Socket,
    {
      conversationUUID,
      restart,
    }: { conversationUUID: string; restart: boolean },
  ) {
    const res = { success: true, data: '' };
    this.chatsService.setCacheObject(
      `${TEXT_ANSWER_GENERATING}:${conversationUUID}`,
      true,
    );
    try {
      await this.chatsService.handleRestartConversation(
        client,
        conversationUUID,
        restart,
      );
    } catch (e) {
      res.success = false;
    }
    return res;
  }
}
