import { Body, Controller, Get, Logger, Param, Post } from '@nestjs/common';
import { LlmsService } from './llms.service';
import { XResponse } from '../../common/response';
import { UpdateLlmsDto } from './dto/update-llms.dto';

@Controller('tools/llms')
export class LlmsController {
  logger = new Logger(LlmsController.name);
  constructor(private readonly llmsService: LlmsService) {}

  @Get()
  async getLLms() {
    return XResponse.Success(await this.llmsService.getLlms());
  }

  @Post('/:uuid')
  async updateLLM(@Param('uuid') uuid: string, @Body() update: UpdateLlmsDto) {
    // this.logger.log(`uuid: ${uuid}`);
    const result = await this.llmsService.updateLLM(uuid, update);
    this.logger.log(`update ${result.affected}`);
    return XResponse.Success('update success');
  }
  @Get('/:uuid')
  async getLLM(@Param('uuid') uuid: string) {
    return XResponse.Success(await this.llmsService.getLLM(uuid));
  }
}
