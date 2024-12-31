import { Module } from '@nestjs/common'
import { UserService } from './user.service'
import { SeedUserService } from './seed-user.service'
import { ToolsService } from '../tools/tools.service';
import { ExperimentsService } from '../experiments/experiments.service';
import { ConversationService } from 'src/chats/conversation.service';


@Module({
  imports: [],
  controllers: [],
  providers: [
    UserService,
    SeedUserService,
    ExperimentsService,
    ToolsService,
    ConversationService,
  ],
  exports: [
    UserService,
  ],
})
export class UserModule {
}