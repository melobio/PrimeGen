import { Injectable, Logger } from '@nestjs/common';
import { DataSource } from 'typeorm';
import { AddUserDto } from './dto/add-user.dto';
import { User } from './entities/user.entity';
import { calculateSHA256Hash } from '../common/utils';
import { ExperimentsService } from '../experiments/experiments.service';
import { CODES, XException } from '../common/exception';

@Injectable()
export class UserService {
  logger = new Logger(UserService.name);
  constructor(
    private readonly dataSource: DataSource,
    private readonly experimentsService: ExperimentsService,
  ) {}

  async findUserByUserName(userName: string) {
    return this.dataSource.manager.findOneBy(User, {
      userName,
    });
  }

  async addUser(addUserDto: AddUserDto) {
    let user = new User();
    user.userName = addUserDto.userName;
    user.email = addUserDto.email;
    user.phone = addUserDto.phone;
    user.password = calculateSHA256Hash(addUserDto.password);

    user = await this.dataSource.manager.save(user);

    return user;
  }

  async registerByUserName(userName: string, password: string) {
    let newUser = new User();
    newUser.userName = userName;
    newUser.password = calculateSHA256Hash(password);
    newUser = await this.dataSource.manager.save(newUser);
    this.experimentsService.addExperiment(newUser, true);
    return newUser;
  }
}
