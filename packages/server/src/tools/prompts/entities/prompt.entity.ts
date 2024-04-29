import {
  Column,
  CreateDateColumn,
  Entity,
  PrimaryGeneratedColumn,
  UpdateDateColumn,
} from 'typeorm';
import { ApiProperty } from '@nestjs/swagger';
import {
  FunctionItem,
  NodeChild,
  NodeField,
  NodeTypes,
  PromptNode,
  PromptType,
} from '@xpcr/common';

@Entity('prompt')
export class PromptEntity implements PromptNode {
  constructor(
    promptType: PromptType,
    name: string,
    shortDesc: string,
    desc: string,
    functions: FunctionItem[],
    fields: NodeField[],
    children: NodeChild[],
    parent: NodeChild,
    version: string,
  ) {
    this.promptType = promptType;
    this.name = name;
    this.shortDesc = shortDesc;
    this.desc = desc;
    this.functions = functions;
    this.fields = fields;
    this.children = children;
    this.parent = parent;
    this.version = version;
  }

  @ApiProperty()
  @PrimaryGeneratedColumn()
  id: number;

  @ApiProperty()
  @Column({ generated: 'uuid' })
  uuid: string;

  @Column()
  name: string;

  type: NodeTypes.Prompts = NodeTypes.Prompts;
  @Column('enum', {
    enum: PromptType,
    default: PromptType.PromptTemplate,
  })
  promptType: PromptType;
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

  @Column('text')
  desc: string;

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
