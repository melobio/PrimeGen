import { Module } from '@nestjs/common';
import { ConversationService } from './conversation.service';
import { ConversationController } from './conversation.controller';
import { ConversationGateway } from './conversation.gateway';

@Module({
  controllers: [ConversationController],
  exports: [ConversationGateway, ConversationService],
  providers: [ConversationService, ConversationGateway],
})
export class ConversationModule {}
