import { IsOptional } from 'class-validator';
import { NodeField } from '@xpcr/common';

export class UpdateLlmsDto {
  @IsOptional()
  fields: NodeField[];
}
