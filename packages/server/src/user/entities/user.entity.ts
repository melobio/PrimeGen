import { Column, CreateDateColumn, Entity, PrimaryGeneratedColumn, UpdateDateColumn } from 'typeorm'
import { ApiProperty } from '@nestjs/swagger'
import { Exclude } from 'class-transformer'

@Entity()
export class User {

  @ApiProperty()
  @PrimaryGeneratedColumn()
  id: number;

  @ApiProperty()
  @Column({ generated: 'uuid' })
  uuid: string;

  @ApiProperty()
  @Column()
  userName: string;

  @ApiProperty()
  @Column({ nullable: true })
  email: string;
  @ApiProperty()
  @Column({ nullable: true })
  phone: string;

  @Column()
  @Exclude()
  password: string;

  @ApiProperty()
  @CreateDateColumn({
    type: 'datetime'
  })
  createTime: Date;

  @ApiProperty()
  @UpdateDateColumn({
    type: 'datetime'
  })
  updateTime: Date;
}