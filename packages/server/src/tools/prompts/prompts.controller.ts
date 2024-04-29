import { Body, Controller, Get, Logger, Param, Post } from '@nestjs/common';
import { PromptsService } from './prompts.service';
import { XResponse } from '../../common/response';
import { PromptEntity } from './entities/prompt.entity';
import { UpdatePromptDto } from './dto/update-prompt.dto';

@Controller('tools/prompts')
export class PromptsController {
  logger = new Logger(PromptsController.name);
  constructor(private readonly promptsService: PromptsService) {}

  @Get()
  async getPrompts() {
    //
    return XResponse.Success<PromptEntity[]>(
      await this.promptsService.getPrompts(),
    );
  }
  @Post('/:uuid')
  async updatePrompt(
    @Param('uuid') uuid: string,
    @Body() update: UpdatePromptDto,
  ) {
    // this.logger.log(`uuid: ${uuid}`);
    const result = await this.promptsService.updatePrompt(uuid, update);
    this.logger.log(`update ${result.affected}`);
    return XResponse.Success('update success');
  }

  @Get('/:uuid')
  getPrompt(@Param('uuid') uuid: string) {
    return XResponse.Success(this.promptsService.getPrompt(uuid));
  }
}
