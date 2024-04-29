import { Inject, Injectable } from '@nestjs/common'
import { CacheService } from '../cache/cache.service'
import { CACHE_MODULE } from '../cache/cache.module'
import { JwtService } from '@nestjs/jwt'
import { UserService } from '../user/user.service'
import { CODES, XException } from '../common/exception'
import { calculateSHA256Hash } from '../common/utils'
import { User } from '../user/entities/user.entity'

@Injectable()
export class AuthService {
  constructor(
    private readonly userService: UserService,
    private readonly jwtService: JwtService,
    @Inject(CACHE_MODULE) private readonly cacheService: CacheService,
    ) {
  }

  async loginByUserName(userName: string, password: string) {
    const user = await this.userService.findUserByUserName(userName);
    if (user?.password !== calculateSHA256Hash(password)) {
      throw new XException(CODES.AUTH.USER_NOT_FOUND_OR_INVALID_PWD, '用户或密码不正确');
    }
    // Generate token
    const token = this.jwtService.sign({
      userName: user.userName,
      uuid: user.uuid,
    })

    await this.cacheService.setToken(token, token);

    return token;
  }
  async verifyToken(token: string) {
    return await this.jwtService.verifyAsync<{ userName: string; uuid: string; }>(token, {  })
  }
  async refreshToken(token: string, user: User) {
    const newToken = this.jwtService.sign({
      userName: user.userName,
      uuid: user.uuid,
    })

    await this.cacheService.setToken(token, newToken);
  }
}