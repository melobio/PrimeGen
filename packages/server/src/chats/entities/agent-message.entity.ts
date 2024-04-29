import {
  Column,
  CreateDateColumn,
  Entity,
  PrimaryGeneratedColumn,
  UpdateDateColumn,
} from 'typeorm';
import { ApiProperty } from '@nestjs/swagger';
import { Role } from './message.entity';

@Entity('agent_messages')
export class AgentMessageEntity {
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
    default: '',
  })
  content: string;
  @Column('text', {
    default: '',
  })
  agentType: string;

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
