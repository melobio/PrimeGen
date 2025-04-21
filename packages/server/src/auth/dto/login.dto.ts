import { ApiProperty } from '@nestjs/swagger'
import { IsEnum, IsNotEmpty, IsOptional, Matches, ValidateIf } from 'class-validator'

export enum LoginType {
  UserName = 'userName',
  Phone = 'phone',
  Email = 'email',
}

export class LoginDto {

  @ApiProperty({ description: '登录类型', enum: LoginType })
  @IsEnum(LoginType)
  loginType: LoginType;

  @ApiProperty({ description: '用户名' })
  @IsNotEmpty()
  @ValidateIf((obj: LoginDto) => obj.loginType === LoginType.UserName)
  userName: string;

  @ApiProperty({ description: '邮箱地址' })
  @Matches(
    /[\w!#$%&'*+/=?^_{|}~-]+(?:.[\w!#$%&'*+/=?^_{|}~-]+)*@(?:[\w](?:[\w-]*[\w])?\.)+[\w](?:[\w-]*[\w])?/,
    {
      message: 'email: 邮箱格式不对'
    }
  )
  @ValidateIf((obj: LoginDto) => obj.loginType === LoginType.Email)
  email: string;

  @ApiProperty({ description: '密码' })
  @IsNotEmpty()
  @ValidateIf((obj: LoginDto) => obj.loginType === LoginType.UserName ||
    obj.loginType === LoginType.Email)
  password: string;

  @ApiProperty({ description: '手机号' })
  @Matches(/^1[34578]\d{9}$/)
  @ValidateIf((obj: LoginDto) => obj.loginType === LoginType.Phone)
  phone: string;

  @ApiProperty({ description: '短信验证码' })
  @Matches(/^\d{6}$/)
  @ValidateIf((obj: LoginDto) => obj.loginType === LoginType.Phone)
  smsCode: string;
}

export class LoginSuccessDto {
  @ApiProperty({ description: 'Token' })
  token: string;
}