import { Module } from '@nestjs/common';
import { ExperimentsService } from './experiments.service';
import { ConversationService } from 'src/chats/conversation.service';
import { ExperimentsController } from './experiments.controller';
import { ToolsModule } from '../tools/tools.module';

@Module({
  controllers: [ExperimentsController],
  providers: [ExperimentsService, ConversationService],
  imports: [ToolsModule],
})
export class ExperimentsModule {}
