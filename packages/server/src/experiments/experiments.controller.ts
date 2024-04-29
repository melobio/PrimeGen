import { Controller, Get, Post, Req } from '@nestjs/common';
import { ExperimentsService } from './experiments.service';
import { XResponse } from '../common/response';
import { Request } from 'express';
import { User } from '../user/entities/user.entity';

@Controller('experiments')
export class ExperimentsController {
  constructor(private readonly experimentsService: ExperimentsService) {}

  @Get()
  async getAllExperiments(@Req() req: Request & { user: User }) {
    return XResponse.Success(
      await this.experimentsService.getAllExperiments(req.user),
    );
  }

  @Post('/addConWithExp')
  async addExperiment(@Req() req: Request & { user: User }) {
    return XResponse.Success(
      await this.experimentsService.addExperiment(req.user),
    );
  }
}
