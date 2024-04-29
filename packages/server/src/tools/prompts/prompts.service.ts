import { Injectable } from '@nestjs/common';
import { DataSource } from 'typeorm';
import { PromptEntity } from './entities/prompt.entity';
import { UpdatePromptDto } from './dto/update-prompt.dto';

@Injectable()
export class PromptsService {
  constructor(private readonly dataSource: DataSource) {}

  async getPrompts() {
    return await this.dataSource.manager.find(PromptEntity);
  }
  async updatePrompt(uuid: string, update: UpdatePromptDto) {
    return await this.dataSource.manager
      .createQueryBuilder()
      .from(PromptEntity, 'prompt')
      .update()
      .set({})
      .where({ uuid })
      .execute();
  }
  async getPrompt(uuid: string) {
    return await this.dataSource.manager.findOne(PromptEntity, {
      where: { uuid },
    });
  }
}
