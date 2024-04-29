import { Injectable, Logger, OnModuleInit } from '@nestjs/common';
import { DataSource } from 'typeorm';
import { User } from './entities/user.entity';
import { calculateSHA256Hash } from '../common/utils';
import * as fs from 'fs';
import * as path from 'path';

@Injectable()
export class SeedUserService implements OnModuleInit {
  logger = new Logger(SeedUserService.name);
  constructor(private readonly dataSource: DataSource) {}

  async create() {
    const users = fs.readFileSync(
      path.join(process.cwd(), 'assets', 'data', 'users.json'),
    );
    const usersArr = JSON.parse(String(users)) as User[];

    for (const user of usersArr) {
      const admin = new User();
      admin.userName = user.userName;
      admin.phone = user.phone;
      admin.email = user.email;
      admin.password = calculateSHA256Hash(user.password);

      const exists = await this.dataSource.manager.exists(User, {
        where: { userName: admin.userName },
      });

      if (!exists) {
        const newAdmin = await this.dataSource.manager.save(admin);
        this.logger.log(`create new Admin: ${JSON.stringify(newAdmin)}`);
      } else {
        this.logger.warn(`${user.userName} exists.`);
      }
    }
  }

  async onModuleInit() {
    this.logger.log(`Start seed admin user.`);
    await this.create();
  }
}
