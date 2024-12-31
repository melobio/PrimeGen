import { ApiProperty } from '@nestjs/swagger'
import { IsEnum, IsNotEmpty, IsOptional, Matches, ValidateIf } from 'class-validator'

export enum LoginType {
  UserName = 'userName',
  Phone = 'phone',
  Email = 'email',
}

export class LoginDto {

  @ApiProperty({ description: 'Login Type', enum: LoginType })
  @IsEnum(LoginType)
  loginType: LoginType;

  @ApiProperty({ description: 'User name' })
  @IsNotEmpty()
  @ValidateIf((obj: LoginDto) => obj.loginType === LoginType.UserName)
  userName: string;

  @ApiProperty({ description: 'Email Address' })
  @Matches(
    /[\w!#$%&'*+/=?^_{|}~-]+(?:.[\w!#$%&'*+/=?^_{|}~-]+)*@(?:[\w](?:[\w-]*[\w])?\.)+[\w](?:[\w-]*[\w])?/,
    {
      message: 'email: invalid email address',
    }
  )
  @ValidateIf((obj: LoginDto) => obj.loginType === LoginType.Email)
  email: string;

  @ApiProperty({ description: 'Password' })
  @IsNotEmpty()
  @ValidateIf((obj: LoginDto) => obj.loginType === LoginType.UserName ||
    obj.loginType === LoginType.Email)
  password: string;

  @ApiProperty({ description: 'Phone' })
  @Matches(/^1[34578]\d{9}$/)
  @ValidateIf((obj: LoginDto) => obj.loginType === LoginType.Phone)
  phone: string;

  @ApiProperty({ description: 'SMS Code' })
  @Matches(/^\d{6}$/)
  @ValidateIf((obj: LoginDto) => obj.loginType === LoginType.Phone)
  smsCode: string;
}

export class RegisterDto {
  @ApiProperty({ description: 'User name' })
  @IsNotEmpty()
  userName: string;

  @ApiProperty({ description: 'Password' })
  @IsNotEmpty()
  password: string;

  @ApiProperty({ description: 'Confirm Password' })
  @IsNotEmpty()
  confirmPassword: string;
}

export class LoginSuccessDto {
  @ApiProperty({ description: 'Token' })
  token: string;
}