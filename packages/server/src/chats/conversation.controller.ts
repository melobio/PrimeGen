import {
  Body,
  Controller,
  Get,
  Logger,
  Param,
  Post,
  Req,
} from '@nestjs/common';
import { ConversationService } from './conversation.service';
import { XResponse } from '../common/response';
import { User } from '../user/entities/user.entity';
import { Request } from 'express';

@Controller('conversations')
export class ConversationController {
  logger = new Logger(ConversationController.name);
  constructor(private readonly chatsService: ConversationService) {}

  @Get()
  async getAllChats(@Req() req: Request & { user: User }) {
    return XResponse.Success(await this.chatsService.getAllChats(req.user));
  }

  @Get(':conversationUUID/messages')
  async getAllMessages(@Param('conversationUUID') conversationUUID: string) {
    return XResponse.Success(
      await this.chatsService.findAllMessages(conversationUUID),
    );
  }

  @Post('/deleteMessageById')
  async deleteMessageById(@Body() messageBo: { id: number }) {
    this.logger.log(
      `ConversationController updateMessageById: ${messageBo.id}`,
    );
    return XResponse.Success(
      await this.chatsService.deleteMessageById(messageBo),
    );
  }

  @Post(':conversationUUID/del')
  async deleteConversationsByUUID(
    @Param('conversationUUID') conversationUUID: string,
  ) {
    return XResponse.Success(
      await this.chatsService.deleteConversationsByUUID(conversationUUID),
    );
  }

  @Post(':conversationUUID/rename')
  async reNameConversationByUUID(
    @Param('conversationUUID') conversationUUID: string,
    @Body() msgInfo: { name: string },
  ) {
    return XResponse.Success(
      await this.chatsService.reNameConversationByUUID(
        conversationUUID,
        msgInfo.name,
      ),
    );
  }
}
