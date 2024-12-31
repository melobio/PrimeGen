import {
  Column,
  CreateDateColumn,
  Entity,
  PrimaryGeneratedColumn,
  UpdateDateColumn,
} from 'typeorm';
import { ApiProperty } from '@nestjs/swagger';
// import { TemplateArgs } from '../../template-args';
import {
  FunctionItem,
  LLMNode,
  LLMType,
  NodeChild,
  NodeField,
} from '@xpcr/common';
import { NodeTypes } from '@xpcr/common';

@Entity('llms')
export class LlmsEntity implements LLMNode {
  constructor(
    llmType: LLMType,
    name: string,
    shortDesc: string,
    desc: string,
    functions: FunctionItem[],
    fields: NodeField[],
    children: NodeChild[],
    parent: NodeChild,
    version: string,
    conversationUUID: string,
  ) {
    this.llmType = llmType;
    this.name = name;
    this.shortDesc = shortDesc;
    this.desc = desc;
    this.functions = functions;
    this.fields = fields;
    this.children = children;
    this.parent = parent;
    this.version = version;
    this.conversationUUID = conversationUUID;

    // this.tempArgs = tempArgs;
  }

  @ApiProperty()
  @PrimaryGeneratedColumn()
  id: number;

  @ApiProperty()
  @Column({ generated: 'uuid' })
  uuid: string;

  @Column()
  name: string;

  type: NodeTypes.LLMs = NodeTypes.LLMs;
  @Column('enum', {
    enum: LLMType,
    default: LLMType.Azure,
  })
  llmType: LLMType;

  @Column('simple-json')
  fields: NodeField[];

  @Column('simple-json')
  functions: FunctionItem[];

  @Column('simple-json')
  children: NodeChild[];

  @Column('simple-json')
  parent: NodeChild;

  @Column('text')
  shortDesc: string;

  @Column()
  conversationUUID: string;

  @Column('text')
  desc: string;

  // @Column('simple-json')
  // tempArgs: { [key: string]: TemplateArgs };

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
