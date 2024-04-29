import { Injectable, Logger, OnModuleInit } from '@nestjs/common';
import { DataSource } from 'typeorm';
import { DevicesEntity } from './entities/devices.entity';
import { UpdateDeviceDto } from './dto/update-device.dto';
import { AgentsService } from '../agents/agents.service';

@Injectable()
export class DevicesService implements OnModuleInit {
  logger = new Logger(DevicesService.name);
  constructor(
    private readonly dataSource: DataSource,
    private readonly agentsService: AgentsService,
  ) {}

  async allDevice() {
    return await this.dataSource.manager.find(DevicesEntity);
  }

  async updateDevice(uuid: string, update: UpdateDeviceDto) {
    return await this.dataSource.manager
      .createQueryBuilder()
      .from(DevicesEntity, 'device')
      .update()
      .set({
        fields: update.fields,
      })
      .where({ uuid })
      .execute();
  }
  async getDevice(uuid: string) {
    return await this.dataSource.manager.findOne(DevicesEntity, {
      where: { uuid },
    });
  }

  // async checkDeviceState(uuid: string) {
  //   const device = await this.dataSource.manager.findOne(DevicesEntity, {
  //     where: { uuid },
  //   });
  //   if (!device) {
  //     throw new XException(CODES.DEVICE.DEVICE_NOT_FOUND, 'Device not found.');
  //   }
  //   this.logger.log(`check device ${device.apiBase} state.`);
  //   let online = false;
  //   try {
  //     const res = (await this.otApi.health(device.apiBase)).data;
  //     this.logger.log(
  //       `check device ${device.apiBase} response: ${JSON.stringify(res)}`,
  //     );
  //     online = Boolean(res);
  //   } catch (e) {
  //     this.logger.error(e, e);
  //   }
  //   // update device state
  //   await this.dataSource.manager
  //     .createQueryBuilder()
  //     .from(DevicesEntity, 'device')
  //     .update()
  //     .set({})
  //     .where({ uuid })
  //     .execute();
  //
  //   return { online };
  // }

  async onModuleInit() {
    // const exists = await this.dataSource.manager.exists(DevicesEntity);
    // if (exists) {
    //   this.logger.warn(`Devices exists.`);
    // } else {
    //   // Seed Devices
    //   const ot2 = await this.dataSource.manager.save(
    //     new DevicesEntity(
    //       'Opentrons OT-2 Robot',
    //       'Description:\n' +
    //         '\n' +
    //         '\n' +
    //         '●OT-2 is a laboratory automation robot designed for high-precision liquid handling and sample analysis.\n' +
    //         '●It is capable of autonomously performing a variety of lab tasks such as pipetting, mixing, centrifugation, and temperature control.\n' +
    //         '●Supports a range of experimental container formats, including test tubes, petri dishes, and microplates.\n' +
    //         '\n' +
    //         'Interface and Communication:\n' +
    //         '\n' +
    //         '\n' +
    //         '●OT-2 communicates with external systems using the TCP/IP protocol.\n' +
    //         '●Supports RESTful API, allowing commands to be sent and responses received via HTTP requests.\n' +
    //         '●Provides command and response formats in JSON, ensuring cross-platform compatibility.\n',
    //       'http://127.0.0.1:31950',
    //       '0.0.1',
    //     ),
    //   );
    //
    //   await this.dataSource.manager.save(
    //     new DevicesEntity(
    //       'MGI AlphaTool Robot',
    //       'Description:\n' +
    //         '\n' +
    //         '\n' +
    //         '●OT-2 is a laboratory automation robot designed for high-precision liquid handling and sample analysis.\n' +
    //         '●It is capable of autonomously performing a variety of lab tasks such as pipetting, mixing, centrifugation, and temperature control.\n' +
    //         '●Supports a range of experimental container formats, including test tubes, petri dishes, and microplates.\n' +
    //         '\n' +
    //         'Interface and Communication:\n' +
    //         '\n' +
    //         '\n' +
    //         '●OT-2 communicates with external systems using the TCP/IP protocol.\n' +
    //         '●Supports RESTful API, allowing commands to be sent and responses received via HTTP requests.\n' +
    //         '●Provides command and response formats in JSON, ensuring cross-platform compatibility.\n',
    //       'http://127.0.0.1:31950',
    //       '0.0.1',
    //     ),
    //   );
    //
    //   setTimeout(async () => {
    //     await this.postInit(ot2);
    //   }, 1000);
    // }
  }

  // private async postInit(device: DevicesEntity) {
  //   // 初始化子节点，OT2作为 AGENT_CODE_EXECUTION_NAME 的子节点
  //   const codeExecution = await this.dataSource.manager.findOne(Agents, {
  //     where: { name: AGENT_CODE_EXECUTION_NAME },
  //   });
  //   if (codeExecution) {
  //     await this.agentsService.updateAgent(codeExecution.uuid, {
  //       // type 需要参考 packages/client/src/pages/main/tools/ToolItem.ts
  //       // 否则会导致渲染失败
  //       children: [],
  //     });
  //     this.logger.log(
  //       `postInit: Agent "${AGENT_CODE_EXECUTION_NAME}" add child ${device.uuid}.`,
  //     );
  //   } else {
  //     this.logger.warn(
  //       `postInit: Agent "${AGENT_CODE_EXECUTION_NAME}" not found.`,
  //     );
  //   }
  // }
}
