import {
  Column,
  CreateDateColumn,
  Entity,
  PrimaryGeneratedColumn,
  UpdateDateColumn,
} from 'typeorm';
import { ApiProperty } from '@nestjs/swagger';
import { NodeChild, ExperimentsState, ExperimentInterface } from '@xpcr/common';

@Entity({ name: 'experiment' })
export class ExperimentEntity implements ExperimentInterface {
  constructor(
    name: string,
    nodeChildren: NodeChild[],
    desc: string,
    creator: string,
    conversationUUID: string,
  ) {
    this.name = name;
    this.nodeChildren = nodeChildren;
    this.desc = desc;
    this.creator = creator;
    this.conversationUUID = conversationUUID;
  }

  @ApiProperty()
  @PrimaryGeneratedColumn()
  id: number;

  @ApiProperty()
  @Column({ generated: 'uuid' })
  uuid: string;

  @Column()
  name: string;

  @Column('simple-json')
  nodeChildren: NodeChild[];

  @Column('text')
  desc: string;

  @Column('uuid')
  creator: string;

  @Column('uuid')
  conversationUUID: string;

  @Column('enum', {
    enum: ExperimentsState,
    default: ExperimentsState.CREATED,
  })
  state: ExperimentsState;

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
