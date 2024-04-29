import {
  MiddlewareConsumer,
  Module,
  NestModule,
  RequestMethod,
} from '@nestjs/common';
import { AppController } from './app.controller';
import { AppService } from './app.service';
import { TypeOrmModule } from '@nestjs/typeorm';
import { User } from './user/entities/user.entity';
import { UserModule } from './user/user.module';
import { CacheModule } from './cache/cache.module';
import { JwtModule } from '@nestjs/jwt';
import { AuthModule } from './auth/auth.module';
import { AllExceptionsFilter } from './common/all-exception-filter';
import { APP_FILTER } from '@nestjs/core';
import { ToolsModule } from './tools/tools.module';
import { Agents } from './tools/agents/entities/agents.entity';
import { AuthMiddleware } from './middleware/auth.middleware';
import { InputEntity } from './tools/input/entities/input.entity';
import { PromptEntity } from './tools/prompts/entities/prompt.entity';
import { LlmsEntity } from './tools/llms/entities/llms.entity';
import { ExperimentsModule } from './experiments/experiments.module';
import { ExperimentEntity } from './experiments/entities/experiment.entity';
import { ConversationModule } from './chats/conversation.module';
import { ConversationEntity } from './chats/entities/conversation.entity';
import { join } from 'path';
import { ServeStaticModule } from '@nestjs/serve-static';
import { MessageEntity } from './chats/entities/message.entity';
import { DevicesEntity } from './tools/devices/entities/devices.entity';
import { ToolEntity } from './tools/tool/entities/tool.entity';
import { FilesModule } from './files/files.module';
import { AgentMessageEntity } from './chats/entities/agent-message.entity';

@Module({
  imports: [
    // 前端
    ServeStaticModule.forRoot({
      rootPath: join(__dirname, '..', '..', 'client/dist'),
      serveRoot: '/',
      exclude: ['/api/(.*)'],
    }),
    JwtModule.register({
      global: true,
      secret: process.env.JWT_SECRET,
      signOptions: { expiresIn: `${process.env.JWT_EXPIRES_IN}s` },
    }),
    CacheModule.register({
      isGlobal: true,
      socket: {
        host: process.env.REDIS_HOST,
        port: +process.env.REDIS_PORT,
      },
      password: process.env.REDIS_PASSWORD,
      database: +process.env.REDIS_DB,
      ttl: +process.env.REDIS_TTL,
    }),
    TypeOrmModule.forRoot({
      type: 'mysql',
      timezone: '+08:00',
      host: process.env.MYSQL_HOST,
      port: +process.env.MYSQL_PORT,
      username: process.env.MYSQL_USERNAME,
      password: process.env.MYSQL_ROOT_PASSWORD,
      database: process.env.MYSQL_DATABASE,
      entities: [
        User,
        Agents,
        InputEntity,
        PromptEntity,
        LlmsEntity,
        ExperimentEntity,
        ConversationEntity,
        MessageEntity,
        AgentMessageEntity,
        DevicesEntity,
        ToolEntity,
      ],
      entityPrefix: process.env.DATABASE_TABLE_PREFIX,
      synchronize: true,
    }),
    AuthModule,
    UserModule,
    ToolsModule,
    ExperimentsModule,
    ConversationModule,
    FilesModule,
  ],
  controllers: [AppController],
  providers: [
    AppService,
    { provide: APP_FILTER, useClass: AllExceptionsFilter },
  ],
})
export class AppModule implements NestModule {
  configure(consumer: MiddlewareConsumer): any {
    consumer
      .apply(AuthMiddleware)
      .exclude(
        { path: '/auth/(.*)', method: RequestMethod.ALL, version: '1' },
        { path: '/files/(.*)', method: RequestMethod.ALL, version: '1' },
      )
      .forRoutes('*');
  }
}
