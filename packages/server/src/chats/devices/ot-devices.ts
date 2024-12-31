import { BaseDevices } from './base-devices';
import { Logger } from '@nestjs/common';
import { OTApi } from './opentrons/OTApi';
import { RunActionType, CheckActionType } from './opentrons/ot-types';
import { JetsonApi } from './jetson/JetsonApi';
import { ConversationService } from '../conversation.service';
import { JETSON_BAD_TIP } from '../conversation.constant';

export class OtDevices extends BaseDevices {
  otApi: OTApi;
  jetsonApi: JetsonApi;
  getLogger(): Logger {
    return new Logger(OtDevices.name);
  }

  constructor(
    readonly apiBase: string,
    readonly jetsonApiBase: string,
    private readonly conversationUUID: string,
    private readonly chatsService: ConversationService,
  ) {
    super();
  }

  async init(): Promise<void> {
    await super.init();
    this.otApi = new OTApi(this.apiBase);
    this.jetsonApi = new JetsonApi(this.jetsonApiBase);
  }

  async createProtocol(protocol: string, protocolKey?: string) {
    const ret = await this.otApi.createProtocol(protocol, protocolKey);
    // console.log('create Response', ret);
    this.logger.debug(
      `create protocol success: ${ret.success} data: ${ret.data}`,
    );
    return ret;
  }
  async createRun(protocolId: string) {
    const ret = await this.otApi.createRun(protocolId);
    // console.log('create Response', ret);
    this.logger.debug(
      `create run success: ${ret.success} data: ${JSON.stringify(ret.data)}`,
    );
    return ret;
  }
  async getRun(runId: string) {
    const ret = await this.otApi.getRun(runId);
    this.logger.debug(
      `get run success: ${ret.success} data: ${JSON.stringify(ret.data)}`,
    );
    return ret;
  }
  async playRun(runId: string) {
    const ret = await this.otApi.runAction(runId, RunActionType.PLAY);
    this.logger.debug(
      `play run success: ${ret.success} data: ${JSON.stringify(ret.data)}`,
    );
    return ret;
  }
  async playPause(runId: string) {
    const ret = await this.otApi.runAction(runId, RunActionType.PAUSE);
    this.logger.debug(
      `play pause success: ${ret.success} data: ${JSON.stringify(ret.data)}`,
    );
    return ret;
  }
  async playStop(runId: string) {
    const ret = await this.otApi.runAction(runId, RunActionType.STOP);
    this.logger.debug(
      `play stop success: ${ret.success} data: ${JSON.stringify(ret.data)}`,
    );
    return ret;
  }
  async wsConnect(runId: string) {
    this.jetsonApi.wsConnectInit();
    this.jetsonApi.ws.on('message', async (jetsonRes: string) => {
      try {
        const text = jetsonRes.toString();
        const res = JSON.parse(text);
        const { event_name, bad_tip, pcr, tip, data, plot } = res;
        if (event_name == 'bad_tip_change' && bad_tip > 0) {
          const playPauseResponse = await this.playPause(runId);
          if (playPauseResponse.success) {
            this.logger.debug(
              `OT2 pause success====> data: ${JSON.stringify(
                playPauseResponse.data,
              )}`,
            );
            // stop jetson and check
            this.jetsonApi.check('stop', 'dropTip', false);
            // ask user to submit command
            this.chatsService.wsClient.emit(JETSON_BAD_TIP, {
              runId,
              conversationUUID: this.conversationUUID,
              bad_tip,
              pcr,
              tip,
              event_name,
              data,
              plot,
            });
          }
        }
      } catch {}
    });
  }
  async getLastRunCommands(runId: string) {
    const ret = await this.otApi.getRunCommands(runId, 1);
    this.logger.debug(
      `get last run commands success: ${ret.success} data: ${JSON.stringify(
        ret.data,
      )}`,
    );
    return ret;
  }
  async getProtocolAnalysis(protocolId: string) {
    const ret = await this.otApi.getProtocolAnalysis(protocolId);
    this.logger.debug(
      `get protocol analysis success: ${ret.success} data: ${JSON.stringify(
        ret.data,
      )}`,
    );
    return ret;
  }
  // -- fault check
  async checkTips(
    action: CheckActionType,
    single = true,
    checkType: 'tips' | 'PCR' | 'start' | 'stop' = 'tips',
  ) {
    const ret = await this.jetsonApi.check<{
      success: boolean;
      data: string;
      plot: [];
    }>(checkType, action, single);
    this.logger.debug(
      `check tips success: ${ret.success} data: ${JSON.stringify(
        JSON.stringify(ret.data),
      )}`,
    );
    try {
      if (ret.data?.plot) {
        this.logger.debug(`check tips ${action} success: data:`);
        this.logger.debug(ret.data.data.substring(0, 30));
        this.logger.debug('check tips success: plot:');
        this.logger.debug(ret.data.plot);
      }
    } catch {}

    return ret;
  }
}
