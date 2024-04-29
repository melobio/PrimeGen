import {
  DeviceType,
  FunctionItem,
  NodeChild,
  NodeField,
  ToolNode,
  ToolType,
} from '@xpcr/common';
import {
  Column,
  CreateDateColumn,
  Entity,
  PrimaryGeneratedColumn,
  UpdateDateColumn,
} from 'typeorm';
import { ApiProperty } from '@nestjs/swagger';
import { NodeTypes } from '@xpcr/common';

@Entity('tool')
export class ToolEntity implements ToolNode {
  constructor(
    toolType: ToolType,
    name: string,
    desc: string,
    functions: FunctionItem[],
    fields: NodeField[],
    children: NodeChild[],
    parent: NodeChild,
    version: string,
  ) {
    this.toolType = toolType;
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

  type: NodeTypes.Tools = NodeTypes.Tools;
  @Column('enum', {
    enum: ToolType,
    default: ToolType.GOOGLE_SEARCH,
  })
  toolType: ToolType;
  @Column('simple-json')
  fields: NodeField[];
  @Column('simple-json')
  functions: FunctionItem[];

  @Column('simple-json')
  children: NodeChild[];

  @Column('simple-json')
  parent: NodeChild;

  @Column()
  name: string;

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
