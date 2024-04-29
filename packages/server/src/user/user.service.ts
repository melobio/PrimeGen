import { Injectable, Logger } from '@nestjs/common'
import { DataSource } from 'typeorm'
import { AddUserDto } from './dto/add-user.dto'
import { User } from './entities/user.entity'
import { calculateSHA256Hash } from '../common/utils'

@Injectable()
export class UserService {
  logger = new Logger(UserService.name)
  constructor(private readonly dataSource: DataSource) {
  }

  async findUserByUserName(userName: string) {
    return this.dataSource.manager.findOneBy(User, {
      userName
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
}