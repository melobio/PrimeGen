export class XException extends Error {
  constructor(readonly code: number, message: string) {
    super(message);
  }
}

export const CODES = {
  SUCCESS: 200,
  AUTH: {
    USER_NOT_FOUND_OR_INVALID_PWD: 10001,
    LOGIN_INVALID: 10002,
    USER_ALREADY_EXIST: 10003,
  },
  DEVICE: {
    DEVICE_NOT_FOUND: 20001,
  },
};
