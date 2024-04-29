import { IsOptional } from 'class-validator';
import { NodeChild } from '@xpcr/common/src/node';

export class UpdateAgentDto {
  @IsOptional()
  children: NodeChild[];
}
