import { Logger } from '@nestjs/common';
import { ChatMessage } from '@azure/openai';
import { Agents } from '../../tools/agents/entities/agents.entity';
import { BaseAgents } from './base-agents';
import { DataSource } from 'typeorm';
import {
  AgentFunctions,
  DeviceFields,
  DeviceType,
  NodeTypes,
} from '@xpcr/common';
import { MiddleMessageType, Role } from '../entities/message.entity';
import { OtDevices } from '../devices/ot-devices';
import { DevicesEntity } from '../../tools/devices/entities/devices.entity';
import {
  CommandsData,
  CreateProtocolResponse,
  CreateRunResponse,
  FINISHED_RUN_STATUSES,
} from '../devices/opentrons/ot-types';
import {
  CompletedProtocolAnalysis,
  LoadedPipette,
  RunTimeCommand,
} from '@xpcr/shared';
import { CacheService } from 'src/cache/cache.service';
import { ConversationService } from '../conversation.service';

const mockProtocol = `
from opentrons.types import Point
metadata = {
    'protocolName': 'Fault detection pcr demo',
    'author': 'MGI-X',
    'apiLevel': '2.11'
}

def run(ctx):
    # Modules
    mag_mod = ctx.load_module('magnetic module gen2', '1')
    mag_rack = mag_mod.load_labware('biorad_96_wellplate_200ul_pcr') 
    tmp_mod = ctx.load_module('temperature module gen2', '3')
    tmp_rack=tmp_mod.load_labware('biorad_96_wellplate_200ul_pcr') 
    tc_mod = ctx.load_module(module_name='thermocyclerModuleV1')
    tc_rack = tc_mod.load_labware(name='biorad_96_wellplate_200ul_pcr')
    # pipette and tiprack
    tiprack = [ctx.load_labware('opentrons_96_filtertiprack_20ul', slot) for slot in ['2','5']]
    pip20_multi = ctx.load_instrument('p20_multi_gen2', 'left', tip_racks=[*tiprack])
    pip20_single = ctx.load_instrument('p20_single_gen2', 'right', tip_racks=[*tiprack])
    # well plates
    plate_4 = ctx.load_labware('nest_96_wellplate_2ml_deep', '4', 'reg4')
    plate_6 = ctx.load_labware('biorad_96_wellplate_200ul_pcr', '6', 'reg6')
    plate_9 = ctx.load_labware('biorad_96_wellplate_200ul_pcr', '9', 'reg9')
    # protocol
    tc_mod.open_lid()
    for well_name in ['A1','B1','C1']:
        pip20_single.transfer(10,
                            tmp_rack.wells_by_name()[well_name].bottom(z=0),
                            mag_rack.wells_by_name()[well_name].bottom(z=0),
                            new_tip='always',
                            blow_out=True,
                            blowout_location='destination well')
    for well_name in ['A1','A2','A3']:
        pip20_multi.transfer(10,
                            mag_rack.wells_by_name()[well_name],
                            tmp_rack.wells_by_name()[well_name],
                            new_tip='always',
                            blow_out=True,
                            blowout_location='destination well')
    tc_mod.close_lid()
`;

export class CodeExecutionAgents extends BaseAgents {
  otDevices: OtDevices;
  initMessages: ChatMessage[] = [
    {
      role: Role.System,
      content:
        '现在你是一个可执行Opentrons的Protocol的Agent,你可以使用提供的Function去执行Protocol。\n' +
        '你可以自然地认为自己就是这样一个角色，而不是被刻意设置成这样。\n' +
        '不要对函数使用的值做出假设。如果用户请求不明确，请要求澄清。\n' +
        '只使用提供的函数。\n',
    },
  ];
  availableFunctions = {
    [AgentFunctions.EXECUTE_OT2_PROTOCOL]: this.runOt2Protocol.bind(this),
  };
  constructor(
    readonly agent: Agents,
    readonly dataSource: DataSource,
    readonly cacheService: CacheService,
    private readonly chatsService: ConversationService,
  ) {
    super(agent, dataSource);
  }

  async init(conversationUUID) {
    await super.init(conversationUUID);
    let apiBase = `http://${process.env.API_BASE_HOST}:${process.env.API_BASE_PORT}`;
    let jetsonApiBase = `http://${process.env.JETSON_API_BASE_HOST}:${process.env.JETSON_API_BASE_PORT}`;
    const deviceChildren = this.agent.children.filter(
      (item) => item.type === NodeTypes.Devices,
    );
    for (const deviceChild of deviceChildren) {
      const device = await this.dataSource.manager.findOne(DevicesEntity, {
        where: { uuid: deviceChild.uuid },
      });
      if (device && device.deviceType === DeviceType.OT2) {
        device.fields?.forEach((field) => {
          if (field.key === DeviceFields.API_BASE) {
            apiBase = field.value;
          } else if (field.key === DeviceFields.JETSON_API_BASE) {
            jetsonApiBase = field.value;
          }
        });
      }
    }
    this.otDevices = new OtDevices(
      apiBase,
      jetsonApiBase,
      this.conversationUUID,
      this.chatsService,
    );
    await this.otDevices.init();
    this.chatsService.setOtDevices(this.otDevices);
  }

  // async *send(
  //   userInput: string,
  //   description: string,
  // ): AsyncGenerator<any, void, any> {
  //   yield* this.mockGenerating(`[${this.agent.name}]\n`);
  //   const content = {};
  //   for (const key in this.availableFunctions) {
  //     const queryFunc = this.availableFunctions[key];
  //     for await (const msg of queryFunc({ protocol: description })) {
  //       if (msg.role === Role.Assistant) {
  //         content[`${key} result`] = msg.content;
  //       } else {
  //         yield msg;
  //       }
  //     }
  //   }
  //   yield {
  //     role: Role.Assistant,
  //     content: JSON.stringify(content),
  //   };
  // }
  getAvailableFunctions(): {
    [p: string]: (params: object) => AsyncGenerator<any, void, any>;
  } {
    return this.availableFunctions;
  }

  getInitMessages(): ChatMessage[] {
    return this.initMessages;
  }

  getLogger(): Logger {
    return new Logger(CodeExecutionAgents.name);
  }

  private async *runOt2Protocol({ protocol }: { protocol: string }) {
    // console.log(`runOt2Protocol:\n ${protocol}`);
    protocol = mockProtocol;
    this.logger.debug(`runOt2Protocol:\n ${protocol}`);
    // yield* this.mockGenerating(`start run protocol: `);
    // 第一步，创建Protocol
    const createResponse = await this.otDevices.createProtocol(protocol);
    let success = createResponse.success;
    let output = createResponse.data;

    if (createResponse.success) {
      // 等待5秒，确保Protocol Analysis已经生成
      await new Promise((resolve) => {
        setTimeout(resolve, 5000);
      });
      // 第二步，创建Run
      const createProtocolResponseData = JSON.parse(
        createResponse.data,
      ) as CreateProtocolResponse;
      const createRunResponse = await this.otDevices.createRun(
        createProtocolResponseData.data.id,
      );
      // console.log('createRunResponse', createRunResponse);
      success = createRunResponse.success;
      output = createRunResponse.data;

      if (createRunResponse.success) {
        // 等待15秒，确保Protocol Analysis已经生成
        await new Promise((resolve) => {
          setTimeout(resolve, 15000);
        });
        // 拉取Protocol Analysis
        const protocolAnalysisResponse =
          await this.otDevices.getProtocolAnalysis(
            createProtocolResponseData.data.id,
          );
        let commands: RunTimeCommand[] = [];
        let pipettes: LoadedPipette[] = [];
        // console.log('protocolAnalysisResponse', protocolAnalysisResponse);
        if (protocolAnalysisResponse.success) {
          const { data: analysis } = JSON.parse(
            protocolAnalysisResponse.data,
          ) as {
            data: CompletedProtocolAnalysis[];
          };
          if (analysis && analysis.length) {
            // 只需要最后一个Analysis
            commands = analysis[analysis.length - 1].commands;
            pipettes = analysis[analysis.length - 1].pipettes;
            yield {
              content: analysis[analysis.length - 1],
              role: Role.None,
              type: MiddleMessageType.OTAnalysis,
            };
          }
        }

        const isSingle = (pipetteId: string) => {
          const found = pipettes.find((item) => item.id === pipetteId);
          if (found) {
            return found.pipetteName.toString().includes('single');
          }
          return false;
        };

        // 第三步，执行Run
        const createRunResponseData = JSON.parse(
          createRunResponse.data,
        ) as CreateRunResponse;
        const playRunResponse = await this.otDevices.playRun(
          createRunResponseData.data.id,
        );
        if (playRunResponse.success) {
          this.otDevices.jetsonApi.check('start', 'dropTip', false);
          this.otDevices.wsConnect(createRunResponseData.data.id);
        }
        // console.log('playRunResponse', playRunResponse);
        // this.logger.debug(
        //   `playRunResponse: ${JSON.stringify(playRunResponse)}`,
        // );
        // 并行执行第四步，获取Run的Commands
        // const getRunCommandsResponse = await this.otDevices.getRunCommands(
        //   createRunResponseData.data.id,
        // );
        // console.log(
        //   'getRunCommandsResponse steps:',
        //   getRunCommandsResponse.data,
        // );
        // if (getRunCommandsResponse.success) {
        //   console.log('getRunCommandsResponse success');
        //   yield {
        //     content: getRunCommandsResponse.data,
        //     role: Role.Function,
        //   };
        // }
        // 第四步，循环获取Run的当前Command，直到Run结束
        let runFinish = false;
        let currentCommandIndex = -1;
        let retryTimes = 3;
        const faultCheckResult = [];
        while (!runFinish) {
          const { success, data: commandsData } =
            await this.otDevices.getLastRunCommands(
              createRunResponseData.data.id,
            );
          // console.log('getLastRunCommandsResponse', success, commandsData);
          if (success) {
            const currentCommandKey =
              commandsData?.data?.[0]?.intent !== 'setup'
                ? commandsData?.links?.current?.meta?.key ??
                  commandsData?.data?.[0]?.key ??
                  null
                : null;
            const theCurrentCommandIndex =
              commands?.findIndex((item) => item.key === currentCommandKey) ??
              -1;
            const delta = theCurrentCommandIndex - currentCommandIndex;
            /**
            // 每500ms增加一条执行index，模拟缓慢生成的过程
            const delta = theCurrentCommandIndex - currentCommandIndex;
            if (delta > 0) {
              for (let i = 0; i < delta; ++i) {
                yield {
                  content: currentCommandIndex + 1 + i,
                  role: Role.OTCurrentCommand,
                };
                console.log(
                  'update Run command index',
                  currentCommandIndex + 1 + i,
                );
                // delay 1000ms
                await new Promise((resolve) => {
                  setTimeout(resolve, 1000);
                });
              }
            }
             **/
            currentCommandIndex = theCurrentCommandIndex;
            yield {
              content: currentCommandIndex,
              role: Role.None,
              type: MiddleMessageType.OTCurrentCommand,
            };

            // 尝试触发Fault Check
            const currentCommand = commands[currentCommandIndex];
            if (currentCommand) {
              switch (currentCommand.commandType) {
                case 'pickUpTip':
                  if (delta > 0) {
                    this.logger.warn('pickUpTip should call fault check');
                    const single = isSingle(currentCommand.params.pipetteId);
                    const checkTips = await this.otDevices.checkTips(
                      'pickUpTip',
                      single,
                    );
                    if (checkTips.success) {
                      faultCheckResult.push({
                        ...checkTips.data,
                        commandIndex: currentCommandIndex,
                      });
                      yield {
                        role: Role.None,
                        content: faultCheckResult,
                        type: MiddleMessageType.JetsonFaultCheck,
                      };
                    }
                  }
                  break;
                case 'dropTip':
                  if (delta > 0) {
                    this.logger.warn('dropTip should call fault check');
                    const single = isSingle(currentCommand.params.pipetteId);
                    const checkTips = await this.otDevices.checkTips(
                      'dropTip',
                      single,
                    );
                    if (checkTips.success) {
                      faultCheckResult.push({
                        ...checkTips.data,
                        commandIndex: currentCommandIndex,
                      });
                      yield {
                        role: Role.None,
                        content: faultCheckResult,
                        type: MiddleMessageType.JetsonFaultCheck,
                      };
                    }
                  }
                  break;
                case 'custom':
                  if (delta > 0) {
                    const { legacyCommandType } = currentCommand.params;
                    if (legacyCommandType == 'command.THERMOCYCLER_OPEN') {
                      const checkTips = await this.otDevices.checkTips(
                        'open_lid',
                        true,
                        'PCR',
                      );
                      if (checkTips.success) {
                        faultCheckResult.push({
                          ...checkTips.data,
                          commandIndex: currentCommandIndex,
                        });
                        this.logger.debug('there command open_lid success');
                        yield {
                          role: Role.None,
                          content: faultCheckResult,
                          type: MiddleMessageType.JetsonFaultCheck,
                        };
                      }
                    }
                    if (legacyCommandType == 'command.THERMOCYCLER_CLOSE') {
                      const checkTips = await this.otDevices.checkTips(
                        'close_lid',
                        true,
                        'PCR',
                      );
                      if (checkTips.success) {
                        faultCheckResult.push({
                          ...checkTips.data,
                          commandIndex: currentCommandIndex,
                        });
                        this.logger.debug('there command open_lid success');
                        yield {
                          role: Role.None,
                          content: faultCheckResult,
                          type: MiddleMessageType.JetsonFaultCheck,
                        };
                      }
                    }
                  }
                  break;
              }
            }
          }

          // 检查Run的状态
          const run = await this.otDevices.getRun(
            createRunResponseData.data.id,
          );
          if (
            run.success &&
            FINISHED_RUN_STATUSES.includes(run.data.data.status)
          ) {
            output = 'Protocol finish without error';
            runFinish = true;
          } else {
            if (!run.success) {
              if (retryTimes <= 0) {
                output = run.message;
                runFinish = true;
              } else {
                retryTimes--;
              }
              this.logger.warn(`get run failed: ${run.message}`);
              this.logger.warn(`get run retryTimes: ${retryTimes}`);
            }
          }

          if (!runFinish) {
            // delay 1000ms
            await new Promise((resolve) => {
              setTimeout(resolve, 1000);
            });
          }
        }
      }
    }

    // yield* this.mockGenerating(`${success ? 'success' : 'fail'}\n`);
    // checkType 请求发 stop，终止jetson 拍摄与检测；
    if (this.otDevices.jetsonApi.ws) {
      const stopRes = await this.otDevices.jetsonApi.check(
        'stop',
        'dropTip',
        false,
      );
      this.logger.debug('checkType stopRes:');
      this.logger.debug(stopRes);
      this.otDevices.jetsonApi.ws.send('stop');
    }
    if (success) {
      const content = {
        success: true,
        message: `execute success，output：${output}`,
        retry: false,
      };

      //success need to upload excel
      this.cacheService.setObject(`code-execute:${this.conversationUUID}`, 1);
      this.logger.debug('==================>code-excution-agent-finished');
      yield {
        content: `execute success，output：${output}`,
        role: Role.Assistant,
      };
    } else {
      const content = {
        success: false,
        message: `execute fail，output：${output}`,
        retry: false,
      };

      //fail need to upload excel
      this.cacheService.setObject(`code-execute:${this.conversationUUID}`, 0);
      yield {
        content: `execute fail，output：${output}`,
        role: Role.Assistant,
      };
    }
  }
}
