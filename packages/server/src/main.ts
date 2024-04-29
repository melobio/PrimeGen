import { NestFactory } from '@nestjs/core';
import { AppModule } from './app.module';
import { Log4jsService } from '@quickts/nestjs-log4js';
import { RequestMethod, ValidationPipe, VersioningType } from '@nestjs/common';
import * as cookieParser from 'cookie-parser';
import { DocumentBuilder, SwaggerModule } from '@nestjs/swagger';

async function bootstrap() {
  const logger = new Log4jsService();
  const app = await NestFactory.create(AppModule);
  app.setGlobalPrefix('api', {
    exclude: [{ path: '/files/(.*)', method: RequestMethod.ALL }],
  });
  app.useGlobalPipes(new ValidationPipe({ transform: true }));
  app.use(cookieParser());
  app.enableVersioning({
    type: VersioningType.URI,
    defaultVersion: '1',
  });
  const config = new DocumentBuilder()
    .setTitle('X PCR')
    .setDescription('X PCR Api Document')
    .setVersion('1.0')
    .addBearerAuth()
    .addTag('auth', '鉴权管理')
    .addTag('agents', 'Agents')
    .build();
  const document = SwaggerModule.createDocument(app, config);
  SwaggerModule.setup('doc', app, document);
  await app.listen(3001);
}
bootstrap().then();
