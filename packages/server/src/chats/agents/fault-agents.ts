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
        'Now you are an Agent checking for errors in OT2 devices, and you can use the provided Function to detect errors in OT2 devices.\n' +
        'You can naturally assume that you are such a role, rather than being deliberately set up to do so.\n' +
        'Do not make assumptions about the values used by the Function. If a user request is unclear, ask for clarification.\n' +
        'Use only the provided functions.',
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
            // check pass
            // yield* this.mockGenerating(`success\n`);
            yield {
              content: 'execute success, output: pickUpTip success',
              role: Role.Assistant,
            };
          } else {
            // check failed
            // yield* this.mockGenerating(`failed\n`);
            yield {
              content: 'execute failed, output: pickUpTip failed',
              role: Role.Assistant,
            };
          }
        } else {
          // might be network error
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
        // mock delaying 2s
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
