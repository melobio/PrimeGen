import { Module } from '@nestjs/common';
import { AgentsModule } from './agents/agents.module';
import { ToolsController } from './tools.controller';
import { InputModule } from './input/input.module';
import { PromptsModule } from './prompts/prompts.module';
import { LlmsModule } from './llms/llms.module';
import { DevicesModule } from './devices/devices.module';
import { ToolsService } from './tools.service';
import { ToolModule } from './tool/tool.module';

@Module({
  controllers: [ToolsController],
  imports: [
    AgentsModule,
    InputModule,
    PromptsModule,
    LlmsModule,
    DevicesModule,
    ToolModule,
  ],
  providers: [ToolsService],
  exports: [ToolsService],
})
export class ToolsModule {}
