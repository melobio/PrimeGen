import { IsOptional } from 'class-validator';
import { NodeField } from '@xpcr/common';

export class UpdateToolDto {
  @IsOptional()
  fields: NodeField[];
}
