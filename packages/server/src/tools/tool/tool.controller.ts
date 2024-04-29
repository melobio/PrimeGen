import { Body, Controller, Get, Logger, Param, Post } from '@nestjs/common';
import { ToolService } from './tool.service';
import { ApiTags } from '@nestjs/swagger';
import { UpdateToolDto } from './dto/update-tool.dto';
import { XResponse } from '../../common/response';
import { ToolEntity } from './entities/tool.entity';
@ApiTags('tools')
@Controller('tools/tools')
export class ToolController {
  logger = new Logger(ToolController.name);
  constructor(private readonly toolService: ToolService) {}

  @Get()
  async allTools() {
    return XResponse.Success<ToolEntity[]>(await this.toolService.allTool());
  }

  @Post('/:uuid')
  async updateTool(@Param('uuid') uuid: string, @Body() update: UpdateToolDto) {
    const result = await this.toolService.updateTool(uuid, update);
    this.logger.log(`update ${result.affected}`);
    return XResponse.Success('update success');
  }
  @Get('/:uuid')
  async getTool(@Param('uuid') uuid: string) {
    return XResponse.Success<ToolEntity>(await this.toolService.getTool(uuid));
  }
}
