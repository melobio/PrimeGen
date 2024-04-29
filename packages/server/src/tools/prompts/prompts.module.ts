import { Module } from '@nestjs/common';
import { PromptsService } from './prompts.service';
import { PromptsController } from './prompts.controller';

@Module({
  controllers: [PromptsController],
  providers: [PromptsService],
})
export class PromptsModule {}
