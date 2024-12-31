export class XResponse<T> {
  constructor(
    readonly code: number,
    readonly message: string,
    readonly success: boolean,
    readonly data?: T,
  ) {
  }
}

export type AdapterRes= {
  code: number,
  msg: string,
  data: string,
  success: boolean,
}

export const CODES = {
  SUCCESS: 200,
  AUTH: {
    USER_NOT_FOUND_OR_INVALID_PWD: 10001,
    LOGIN_INVALID: 10002,
  },
  DEVICE: {
    DEVICE_NOT_FOUND: 20001,
  },
};