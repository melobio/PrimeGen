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
      this.logger.warn(`No token in header.`, '');
      throw new XException(CODES.AUTH.LOGIN_INVALID, 'please login again!');
    }
    const cacheToken = await this.cacheService.getToken<string>(token);
    // this.logger.debug(`cacheToken: ${cacheToken}`);
    if (!cacheToken) {
      // throw new UnauthorizedException(`please login again!');`);
      this.logger.warn(`Token #${token} not in cache`, '');
      throw new XException(CODES.AUTH.LOGIN_INVALID, 'please login again!');
    }

    const { userName } = await this.authService.verifyToken(cacheToken);
    const user = await this.userService.findUserByUserName(userName);
    if (!user) {
      this.logger.warn(`User ${userName} not found.`, '');
      throw new XException(CODES.AUTH.LOGIN_INVALID, 'please login again!');
    }

    // refresh cache toke
    await this.authService.refreshToken(token, user);

    // store user in request
    req.user = user;

    next();
  }

  private extractTokenFromHeader(request: Request): string | undefined {
    const [type, token] = request.headers.authorization?.split(' ') ?? [];
    return type === 'Bearer' ? token : undefined;
  }
}
