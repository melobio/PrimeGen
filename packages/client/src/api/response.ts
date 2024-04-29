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