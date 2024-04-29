import { Module } from '@nestjs/common';
import { AgentsController } from './agents.controller';
import { AgentsService } from './agents.service';
import { SeedAgentsService } from './seed-agents.service';

@Module({
  controllers: [AgentsController],
  providers: [AgentsService, SeedAgentsService],
  exports: [AgentsService],
})
export class AgentsModule {}
