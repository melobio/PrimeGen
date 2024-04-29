import { Body, Controller, Post, Res } from '@nestjs/common'
import { ApiTags } from '@nestjs/swagger'
import { LoginDto, LoginSuccessDto, LoginType } from './dto/login.dto'
import { AuthService } from './auth.service'
import { XResponse } from '../common/response'
import { Response } from 'express'

@ApiTags('auth')
@Controller('auth')
export class AuthController {
  constructor(private readonly authService: AuthService) {
  }

  @Post('login')
  async login(
    @Body() loginDto: LoginDto,
    @Res({ passthrough: true }) response: Response,
  ) {
    const loginSuccessDto = new LoginSuccessDto()
    switch (loginDto.loginType) {
      case LoginType.UserName:
        loginSuccessDto.token = await this.authService
          .loginByUserName(loginDto.userName, loginDto.password);
        response.cookie('x-token', loginSuccessDto.token, {
          httpOnly: false,
        });
        return XResponse.Success(loginSuccessDto);
      default:
        return XResponse.Fail(0, `未实现登录方式：${loginDto.loginType}`)
    }
  }
}