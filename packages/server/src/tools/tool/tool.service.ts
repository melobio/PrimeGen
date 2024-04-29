import { Injectable, Logger } from '@nestjs/common';
import { DataSource } from 'typeorm';
import { ToolEntity } from './entities/tool.entity';
import { UpdateToolDto } from './dto/update-tool.dto';

@Injectable()
export class ToolService {
  logger = new Logger(ToolService.name);

  constructor(private readonly dataSource: DataSource) {}

  async allTool() {
    return await this.dataSource.manager.find(ToolEntity);
  }

  async updateTool(uuid: string, update: UpdateToolDto) {
    return await this.dataSource.manager
      .createQueryBuilder()
      .from(ToolEntity, 'tool')
      .update()
      .set({
        fields: update.fields,
      })
      .where({ uuid })
      .execute();
  }

  async getTool(uuid: string) {
    return await this.dataSource.manager.findOne(ToolEntity, {
      where: { uuid },
    });
  }
}
