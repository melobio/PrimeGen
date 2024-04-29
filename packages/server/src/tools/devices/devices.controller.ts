import { Body, Controller, Get, Logger, Param, Post } from '@nestjs/common';
import { DevicesService } from './devices.service';
import { ApiTags } from '@nestjs/swagger';
import { DevicesEntity } from './entities/devices.entity';
import { XResponse } from '../../common/response';
import { UpdateDeviceDto } from './dto/update-device.dto';

@ApiTags('devices')
@Controller('tools/devices')
export class DevicesController {
  logger = new Logger(DevicesController.name);
  constructor(private readonly devicesService: DevicesService) {}

  @Get()
  async allDevices() {
    return XResponse.Success<DevicesEntity[]>(
      await this.devicesService.allDevice(),
    );
  }
  @Post('/:uuid')
  async updateDevice(
    @Param('uuid') uuid: string,
    @Body() update: UpdateDeviceDto,
  ) {
    // this.logger.log(`uuid: ${uuid}`);
    const result = await this.devicesService.updateDevice(uuid, update);
    this.logger.log(`update ${result.affected}`);
    return XResponse.Success('update success');
  }
  @Get('/:uuid')
  async getDevice(@Param('uuid') uuid: string) {
    return XResponse.Success(await this.devicesService.getDevice(uuid));
  }

  // @Get('/:uuid/state')
  // async checkDeviceState(@Param('uuid') uuid: string) {
  //   const result = await this.devicesService.checkDeviceState(uuid);
  //   return XResponse.Success(result);
  // }
}
