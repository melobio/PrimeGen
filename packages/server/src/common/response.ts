import { CODES } from './exception'

export class XResponse<T> {
  constructor(
    private code: number,
    private message: string,
    private success: boolean,
    private data?: T,
  ) {
  }

  static Success<T>(data: T): XResponse<T> {
    return new XResponse<T>(
      CODES.SUCCESS,
      'Success',
      true,
      data
    )
  }

  static Fail(code: number, message: string): XResponse<any> {
    return new XResponse<any>(
      code,
      message,
      false,
      null,
    )
  }
}