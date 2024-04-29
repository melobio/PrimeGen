import { IsOptional } from 'class-validator';
import { NodeField } from '@xpcr/common';

export class UpdateDeviceDto {
  @IsOptional()
  fields: NodeField[];
}
