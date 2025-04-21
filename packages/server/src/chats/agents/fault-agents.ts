import { Logger } from '@nestjs/common';
import { Agents } from '../../tools/agents/entities/agents.entity';
import { ChatMessage } from '@azure/openai';
import { BaseAgents } from './base-agents';
import {
  AgentFunctions,
  DeviceFields,
  DeviceType,
  NodeTypes,
} from '@xpcr/common';
import { DataSource } from 'typeorm';
import { MiddleMessageType, Role } from '../entities/message.entity';
import { OtDevices } from '../devices/ot-devices';
import { DevicesEntity } from '../../tools/devices/entities/devices.entity';
import { ConversationService } from '../conversation.service';

export class FaultAgents extends BaseAgents {
  otDevices: OtDevices;
  initMessages: ChatMessage[] = [
    {
      role: 'system',
      content:
        '现在你是一个检查OT2设备是否出错的Agent,你可以使用提供的Function去检测OT2设备的错误。\n' +
        '你可以自然地认为自己就是这样一个角色，而不是被刻意设置成这样。\n' +
        '不要对函数使用的值做出假设。如果用户请求不明确，请要求澄清。\n' +
        '只使用提供的函数。',
    },
  ];
  availableFunctions = {
    [AgentFunctions.CHECK_OT2_STATE]: this.checkOT2fault.bind(this),
  };
  constructor(
    readonly agent: Agents,
    readonly dataSource: DataSource,
    private readonly chatsService: ConversationService,
  ) {
    super(agent, dataSource);
  }

  async init(conversationUUID) {
    await super.init(conversationUUID);
    let apiBase = 'http://127.0.0.1:31950';
    let jetsonApiBase = 'http://127.0.0.1';
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
  }

  getAvailableFunctions(): {
    [p: string]: (params: object) => AsyncGenerator<any, void, any>;
  } {
    return this.availableFunctions;
  }

  getLogger(): Logger {
    return new Logger(FaultAgents.name);
  }
  getInitMessages(): ChatMessage[] {
    return this.initMessages;
  }
  private async *checkOT2fault({
    checkType = 'all',
  }: {
    checkType: 'all' | 'liquid' | 'tips';
  }) {
    this.logger.debug(`checkOT2fault: ${checkType}`);
    // yield* this.mockGenerating(`start check ${checkType}: `);

    switch (checkType) {
      case 'tips':
        const checkTips = await this.otDevices.checkTips('pickUpTip');
        if (checkTips.success) {
          yield {
            role: Role.None,
            type: MiddleMessageType.JetsonFaultCheck,
            content: [checkTips.data],
          };
          if (checkTips.data.success) {
            // 检查通过
            // yield* this.mockGenerating(`success\n`);
            yield {
              content: 'execute success, output: pickUpTip success',
              role: Role.Assistant,
            };
          } else {
            // 检查不通过
            // yield* this.mockGenerating(`failed\n`);
            yield {
              content: 'execute failed, output: pickUpTip failed',
              role: Role.Assistant,
            };
          }
        } else {
          // 可能是网络失败
          // yield* this.mockGenerating(`failed\n`);
          yield {
            content:
              checkTips.message ?? 'execute failed, please check network.',
            role: Role.Assistant,
          };
        }
        break;
      // case 'liquid':
      //   break;
      // case 'all':
      //   break;
      default:
        // 模拟延时2s
        await new Promise((resolve) => setTimeout(resolve, 2000));
        // yield* this.mockGenerating(`success\n`);
        yield {
          content: 'execute success，output：success',
          role: Role.Assistant,
        };
        break;
    }
  }
}
