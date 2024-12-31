import {
  Column,
  CreateDateColumn,
  Entity,
  PrimaryGeneratedColumn,
  UpdateDateColumn,
} from 'typeorm';
import { ApiProperty } from '@nestjs/swagger';
import {
  DeviceNode,
  DeviceType,
  FunctionItem,
  NodeChild,
  NodeField,
  NodeTypes,
} from '@xpcr/common';

@Entity('devices')
export class DevicesEntity implements DeviceNode {
  constructor(
    deviceType: DeviceType,
    name: string,
    desc: string,
    functions: FunctionItem[],
    fields: NodeField[],
    children: NodeChild[],
    parent: NodeChild,
    version: string,
    conversationUUID: string,
  ) {
    this.deviceType = deviceType;
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

  type: NodeTypes.Devices = NodeTypes.Devices;
  @Column('enum', {
    enum: DeviceType,
    default: DeviceType.OT2,
  })
  deviceType: DeviceType;
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

  @Column()
  conversationUUID: string;

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
