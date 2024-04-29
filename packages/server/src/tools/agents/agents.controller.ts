import { Body, Controller, Get, Logger, Param, Post } from '@nestjs/common';
import { AgentsService } from './agents.service';
import { ApiTags } from '@nestjs/swagger';
import { XResponse } from '../../common/response';
import { Agents } from './entities/agents.entity';
import { UpdateAgentDto } from './dto/update-agent.dto';

@ApiTags('agents')
@Controller('tools/agents')
export class AgentsController {
  logger = new Logger(AgentsController.name);
  constructor(private readonly agentsService: AgentsService) {}

  @Get()
  async allAgents() {
    return XResponse.Success<Agents[]>(await this.agentsService.allAgent());
  }

  @Post('/:uuid')
  async updateAgent(
    @Param('uuid') uuid: string,
    @Body() update: UpdateAgentDto,
  ) {
    // this.logger.log(`uuid: ${uuid}`);
    const result = await this.agentsService.updateAgent(uuid, update);
    this.logger.log(`update ${result.affected}`);
    return XResponse.Success('update success');
  }
  @Get('/:uuid')
  async getAgent(@Param('uuid') uuid: string) {
    return XResponse.Success(await this.agentsService.getAgent(uuid));
  }
}
