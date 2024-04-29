import { ArgumentsHost, BadRequestException, Catch, ExceptionFilter, Logger } from '@nestjs/common'
import { XException } from './exception'
import { Request, Response } from 'express'
import { XResponse } from './response'

@Catch()
export class AllExceptionsFilter implements ExceptionFilter {
  logger = new Logger(AllExceptionsFilter.name);

  catch(exception: any, host: ArgumentsHost): any {
    const ctx = host.switchToHttp();
    const response = ctx.getResponse() as Response;
    const request = ctx.getRequest() as Request;

    this.logger.error(`body ${JSON.stringify(request.body)}`, exception.stack)

    if (exception instanceof XException) {
      const exp = exception as XException;
      response.json(
        XResponse.Fail(exp.code, exp.message)
      )
    } else if (exception instanceof BadRequestException) {
      const exp = exception as BadRequestException;
      const resp = exp.getResponse()
      response
        .json(XResponse.Fail(0, JSON.stringify(resp)));
    } else {
      response
        .json(XResponse.Fail(0, exception.message || 'Unknown Error'));
    }
  }
}