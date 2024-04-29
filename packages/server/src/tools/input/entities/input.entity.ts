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
  InputNode,
  InputType,
  NodeChild,
  NodeField,
  NodeTypes,
} from '@xpcr/common';

@Entity('input')
export class InputEntity implements InputNode {
  constructor(
    inputType: InputType,
    name: string,
    desc: string,
    functions: FunctionItem[],
    fields: NodeField[],
    children: NodeChild[],
    parent: NodeChild,
    version: string,
  ) {
    this.inputType = inputType;
    this.name = name;
    this.desc = desc;
    this.shortDesc = desc;
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

  type: NodeTypes.Input = NodeTypes.Input;
  @Column('enum', {
    enum: InputType,
    default: InputType.InputNode,
  })
  inputType: InputType;
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
