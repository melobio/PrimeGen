import { Controller, Get, Param } from '@nestjs/common';
import { InputService } from './input.service';
import { ApiTags } from '@nestjs/swagger';
import { XResponse } from '../../common/response';
import { InputEntity } from './entities/input.entity';

@ApiTags('inputs')
@Controller('tools/inputs')
export class InputController {
  constructor(private readonly inputService: InputService) {}

  @Get()
  async allInputs() {
    return XResponse.Success<InputEntity[]>(await this.inputService.allInput());
  }

  @Get('/:uuid')
  getInput(@Param('uuid') uuid: string) {
    return XResponse.Success(this.inputService.getInput(uuid));
  }
}
