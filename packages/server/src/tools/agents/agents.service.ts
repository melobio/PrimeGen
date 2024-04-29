import { Injectable } from '@nestjs/common';
import { DataSource } from 'typeorm';
import { Agents } from './entities/agents.entity';
import { UpdateAgentDto } from './dto/update-agent.dto';

@Injectable()
export class AgentsService {
  constructor(private readonly dataSource: DataSource) {}

  async allAgent() {
    return await this.dataSource.manager.find(Agents);
  }

  async updateAgent(uuid: string, update: UpdateAgentDto) {
    return await this.dataSource.manager
      .createQueryBuilder()
      .from(Agents, 'agent')
      .update()
      .set({
        children: update.children,
      })
      .where({ uuid })
      .execute();
  }
  async getAgent(uuid: string) {
    return await this.dataSource.manager.findOne(Agents, {
      where: { uuid },
    });
  }
}
