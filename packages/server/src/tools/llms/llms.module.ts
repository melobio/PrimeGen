import { Module } from '@nestjs/common';
import { LlmsService } from './llms.service';
import { LlmsController } from './llms.controller';

@Module({
  controllers: [LlmsController],
  providers: [LlmsService],
})
export class LlmsModule {}
