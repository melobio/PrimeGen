import { ApiProperty } from '@nestjs/swagger';
import {
  Column,
  CreateDateColumn,
  Entity,
  PrimaryGeneratedColumn,
  UpdateDateColumn,
} from 'typeorm';
import { PlanStep } from '../agents/ext/planner';

@Entity('conversation')
export class ConversationEntity {
  constructor(name: string, desc: string, creator: string) {
    this.name = name;
    this.desc = desc;
    this.creator = creator;
  }
  @ApiProperty()
  @PrimaryGeneratedColumn()
  id: number;

  @ApiProperty()
  @Column({ generated: 'uuid' })
  uuid: string;

  @Column('uuid')
  creator: string;

  @Column()
  name: string;

  @Column('text')
  desc: string;

  @Column('varchar', {
    default: '',
  })
  currentStep: string;

  @Column('json', { nullable: true })
  stepsResult: {
    [stepName: string]: string;
  };

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
