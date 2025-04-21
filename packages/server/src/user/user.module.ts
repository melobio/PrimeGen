import { Module } from '@nestjs/common'
import { UserService } from './user.service'
import { SeedUserService } from './seed-user.service'

@Module({
  imports: [],
  controllers: [],
  providers: [
    UserService,
    SeedUserService,
  ],
  exports: [
    UserService,
  ],
})
export class UserModule {
}