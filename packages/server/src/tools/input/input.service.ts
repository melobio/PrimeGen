import { Injectable } from '@nestjs/common';
import { DataSource } from 'typeorm';
import { InputEntity } from './entities/input.entity';

@Injectable()
export class InputService {
  constructor(private readonly dataSource: DataSource) {}
  async allInput() {
    return this.dataSource.manager.find(InputEntity);
  }
  async getInput(uuid: string) {
    return this.dataSource.manager.findOne(InputEntity, {
      where: { uuid },
    });
  }
}
