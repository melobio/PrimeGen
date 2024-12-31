import {
  Column,
  CreateDateColumn,
  Entity,
  PrimaryGeneratedColumn,
  UpdateDateColumn,
} from 'typeorm';
import { NodeField } from '@xpcr/common';
import { ApiProperty } from '@nestjs/swagger';
import type { AgentNode } from '@xpcr/common/src/agent';
import { AgentType } from '@xpcr/common';
import type { FunctionItem } from '@xpcr/common/src/functions';
import { NodeChild, NodeTypes } from '@xpcr/common';

@Entity()
export class Agents implements AgentNode {
  constructor(
    agentType: AgentType,
    name: string,
    desc: string,
    functions: FunctionItem[],
    fields: NodeField[],
    children: NodeChild[],
    parent: NodeChild,
    version: string,
    conversationUUID: string,
  ) {
    this.agentType = agentType;
    this.name = name;
    this.desc = desc;
    this.shortDesc = desc;
    this.functions = functions;
    this.fields = fields;
    this.children = children;
    this.parent = parent;
    this.version = version;
    this.conversationUUID = conversationUUID;
  }

  @ApiProperty()
  @PrimaryGeneratedColumn()
  id: number;

  @ApiProperty()
  @Column({ generated: 'uuid' })
  uuid: string;

  type: NodeTypes.Agents = NodeTypes.Agents;

  @ApiProperty()
  @Column('enum', {
    enum: AgentType,
    default: AgentType.FAULT,
  })
  agentType: AgentType;

  @Column()
  name: string;

  @Column()
  conversationUUID: string;

  @Column('text')
  shortDesc: string;

  @Column('text')
  desc: string;

  @Column('simple-json')
  fields: NodeField[];
  @Column('simple-json')
  functions: FunctionItem[];

  @Column('simple-json')
  children: NodeChild[];

  @Column('simple-json')
  parent: NodeChild;

  @Column()
  version: string;

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
