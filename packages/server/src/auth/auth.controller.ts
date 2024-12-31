import { Body, Controller, Post, Res } from '@nestjs/common'
import { ApiTags } from '@nestjs/swagger'
import { LoginDto, RegisterDto, LoginSuccessDto, LoginType } from './dto/login.dto'
import { AuthService } from './auth.service'
import { XResponse } from '../common/response'
import { Response } from 'express'

@ApiTags('auth')
@Controller('auth')
export class AuthController {
  constructor(private readonly authService: AuthService) {}

  @Post('login')
  async login(
    @Body() loginDto: LoginDto,
    @Res({ passthrough: true }) response: Response,
  ) {
    const loginSuccessDto = new LoginSuccessDto();
    switch (loginDto.loginType) {
      case LoginType.UserName:
        loginSuccessDto.token = await this.authService.loginByUserName(
          loginDto.userName,
          loginDto.password,
        );
        response.cookie('x-token', loginSuccessDto.token, {
          httpOnly: false,
        });
        return XResponse.Success(loginSuccessDto);
      default:
        return XResponse.Fail(
          0,
          `unsupported login type:${loginDto.loginType}`,
        );
    }
  }

  @Post('register')
  async register(
    @Body() register_body: RegisterDto,
    @Res({ passthrough: true }) response: Response,
  ) {
    register_body.userName = register_body.userName.trim();
    if (register_body.userName.length === 0) {
      return XResponse.Fail(0, 'Username cannot be empty');
    }
    if (register_body.password.length === 0) {
      return XResponse.Fail(0, 'Password cannot be empty');
    }
    // check password length
    if (register_body.password.length < 6) {
      return XResponse.Fail(0, 'Password length is less than 6');
    }
    // require at least one lowercase letter, one uppercase letter, one digit
    const reg = /^(?=.*[a-z])(?=.*[A-Z])(?=.*\d)[\s\S]{6,}$/;
    if (!reg.test(register_body.password)) {
      return XResponse.Fail(
        0,
        `Password must contain: \n
1. At least one lowercase letter \n
2. At least one uppercase letter \n
3. At least one digit \n
4. At least 6 characters \n
5. Can contain any characters (spaces, special chars, etc.) 
        `,
      );
    }
    // check password and confirmPassword
    if (register_body.password !== register_body.confirmPassword) {
      return XResponse.Fail(0, 'Two password entries are inconsistent');
    }
    // check username if exists
    const user = await this.authService.getUser(register_body.userName);
    if (user) {
      return XResponse.Fail(0, 'Username already exists');
    }
    // register
    const newUser = await this.authService.registerByUserName(
      register_body.userName,
      register_body.password,
    );
    return XResponse.Success({
      userName: newUser.userName,
    });
  }
}
