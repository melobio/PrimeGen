import {
  Column,
  CreateDateColumn,
  Entity,
  PrimaryGeneratedColumn,
  UpdateDateColumn,
} from 'typeorm';
import { ApiProperty } from '@nestjs/swagger';
import { RunCommandSummary } from '../devices/opentrons/ot-types';
import { CompletedProtocolAnalysis } from '@xpcr/shared';
import { NCBIFunctionType, PrimerDesignFunctionType } from '../agents/ext/ncbi-search';

export enum Role {
  User = 'user',
  System = 'system',
  Assistant = 'assistant',
  Function = 'function',
  // Agent = 'agent',

  OTAnalysis = 'ot_analysis',
  OTCurrentCommand = 'current_command',
  JetsonFaultCheck = 'jetson_fault_check',

  None = 'none', // 中间消息
}
export enum MiddleMessageType {
  OTAnalysis = 'ot_analysis',
  OTCurrentCommand = 'current_command',
  JetsonFaultCheck = 'jetson_fault_check',
  SequenceSearchOption = 'sequence_search_option',
}

@Entity('messages')
export class MessageEntity {
  @ApiProperty()
  @PrimaryGeneratedColumn()
  id: number;

  @ApiProperty()
  @Column({ generated: 'uuid' })
  uuid: string;

  @Column('uuid')
  conversationUUID: string;

  @Column('enum', {
    enum: Role,
    default: Role.User,
  })
  role: Role;

  @Column('text', {
    // default: '',
  })
  header: string;

  @Column('text', {
    // default: '',
  })
  agentType: string;

  @Column('text')
  content: string;

  @Column('json', {
    nullable: true,
  })
  protocolAnalysis?: CompletedProtocolAnalysis;

  @Column('json', {
    nullable: true,
  })
  optionInfo?: {
    options?: Record<string, any>[];
    cancer_list?: Record<string, any>;
    data_list?: Record<string, any>[];
    selected_options?: string[];
    stage: number;
    state: string;
    search_type?: NCBIFunctionType;
    species_name?: string;
    submitted?: boolean;
    strain_path?: string[];
    non_target_path?: string[];
    target_gene_path?: string[];
    upload_file?: boolean;
    choice_type?: string;
    species_identification_dict?: Record<string, string>;
    protein_mutation_dict?: Record<string, string>;
    primer_design_prompt?: string;
    primer_design_dict?: Record<string, string>;
    primer_type?: PrimerDesignFunctionType | string;
    operations?: {
      title: string;
      key: string;
      options: string[];
      value: any[];
      type: string[];
    }[];
  };

  @Column('simple-json', { nullable: true })
  uploadFiles?: string[];

  @Column({
    nullable: true,
    default: -1,
  })
  currentCommandIndex?: number;

  @Column('json', {
    nullable: true,
  })
  faultCheckResult?: {
    success: boolean;
    data: string;
    plot: number[];
    commandIndex?: number;
  }[];

  @ApiProperty()
  @CreateDateColumn({
    type: 'datetime',
  })
  createTime: Date;

  @ApiProperty()
  @UpdateDateColumn({
    type: 'datetime',
  })
  updateTime: Date;
}
