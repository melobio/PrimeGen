import {
  Inject,
  Injectable,
  Logger,
  NestMiddleware,
  UnauthorizedException,
} from '@nestjs/common';
import { NextFunction, Request, Response } from 'express';
import { UserService } from '../user/user.service';
import { AuthService } from '../auth/auth.service';
import { CACHE_MODULE } from '../cache/cache.module';
import { CacheService } from '../cache/cache.service';
import { JwtService } from '@nestjs/jwt';
import { CODES, XException } from '../common/exception';

@Injectable()
export class AuthMiddleware implements NestMiddleware {
  private readonly logger = new Logger(AuthMiddleware.name);
  constructor(
    private authService: AuthService,
    private readonly jwtService: JwtService,
    private readonly userService: UserService,
    @Inject(CACHE_MODULE) private cacheService: CacheService,
  ) {}
  async use(req: Request, res: Response, next: NextFunction) {
    this.logger.log(`req: ${req.baseUrl}`); //
    const token = this.extractTokenFromHeader(req);
    if (!token) {
      // throw new UnauthorizedException('token不存在');
      this.logger.warn(`No token in header.`, '');
      throw new XException(CODES.AUTH.LOGIN_INVALID, '登录失效，请重新登录');
    }
    const cacheToken = await this.cacheService.getToken<string>(token);
    // this.logger.debug(`cacheToken: ${cacheToken}`);
    if (!cacheToken) {
      // throw new UnauthorizedException(`登录失效，请重新登录`);
      this.logger.warn(`Token #${token} not in cache`, '');
      throw new XException(CODES.AUTH.LOGIN_INVALID, '登录失效，请重新登录');
    }

    const { userName } = await this.authService.verifyToken(cacheToken);
    const user = await this.userService.findUserByUserName(userName);
    if (!user) {
      // throw '用户不存在';
      this.logger.warn(`User ${userName} not found.`, '');
      throw new XException(CODES.AUTH.LOGIN_INVALID, '登录失效，请重新登录');
    }

    // 刷新缓存中的Token，续期
    await this.authService.refreshToken(token, user);

    // 方便后续可以直接获取到登录的用户数据
    req.user = user;

    next();
  }

  private extractTokenFromHeader(request: Request): string | undefined {
    const [type, token] = request.headers.authorization?.split(' ') ?? [];
    return type === 'Bearer' ? token : undefined;
  }
}
