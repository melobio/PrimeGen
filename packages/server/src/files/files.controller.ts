import {
  Controller,
  Get,
  Logger,
  Param,
  Req,
  Res,
  StreamableFile,
} from '@nestjs/common';
import { FilesService } from './files.service';
import { join } from 'path';
import { createReadStream } from 'fs';
import type { Response, Request } from 'express';
import * as fs from 'fs';

@Controller('files')
export class FilesController {
  rootDir = process.env.XPCR_FILES;
  logger = new Logger(FilesController.name);
  constructor(private readonly filesService: FilesService) {}

  @Get('*')
  getFile(
    @Res({ passthrough: true }) res: Response,
    @Req() request: Request,
  ): StreamableFile {
    const filePath = request.params[0].trim();
    const fullPath = join(this.rootDir, filePath);
    this.logger.debug(`get File fullPath:${fullPath}`);
    if (!fs.existsSync(fullPath)) {
      res.status(500);
      throw new Error(`file not found: ${fullPath}`);
    }
    const file = createReadStream(fullPath);
    // res.set({
    //   'Content-Type': 'application/json',
    //   'Content-Disposition': 'attachment; filename="package.json"',
    // });
    return new StreamableFile(file);
  }
}
