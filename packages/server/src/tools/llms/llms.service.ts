import { Injectable, Logger, OnModuleInit } from '@nestjs/common';
import { DataSource } from 'typeorm';
import { LlmsEntity } from './entities/llms.entity';
import { TemplateArgType } from '../template-args';
import { UpdateLlmsDto } from './dto/update-llms.dto';

export const OpenAIModels = [
  'gpt-3.5-turbo-0613',
  'gpt-4-0613',
  'gpt-4-32k-0613',
  'gpt-4',
  'gpt-4-32',
];

@Injectable()
export class LlmsService implements OnModuleInit {
  logger = new Logger(LlmsService.name);
  constructor(private readonly dataSource: DataSource) {}

  async getLlms() {
    return await this.dataSource.manager.find(LlmsEntity);
  }

  async onModuleInit() {
    // const exists = await this.dataSource.manager.exists(LlmsEntity);
    // if (exists) {
    //   this.logger.warn(`LLms exists`);
    // } else {
    //   await this.dataSource.manager.save(
    //     new LlmsEntity(
    //       'ChatOpenAI',
    //       'Wrapper around OpenAI Chat Large language models.',
    //       'Wrapper around OpenAI Chat Large language models.',
    //       {
    //         model_name: {
    //           label: 'Model Name',
    //           key: 'model_name',
    //           value: 'gpt-3.5-turbo-0613',
    //           type: TemplateArgType.Select,
    //           options: OpenAIModels,
    //           hint: '',
    //           desc: '',
    //         },
    //         api_base: {
    //           label: 'OpenAI API Base',
    //           key: 'api_base',
    //           value: '',
    //           type: TemplateArgType.TextField,
    //           hint: 'Input please',
    //           desc: '',
    //         },
    //         api_key: {
    //           label: 'OpenAI API Key',
    //           key: 'api_key',
    //           value: '',
    //           type: TemplateArgType.Secret,
    //           hint: 'Input please',
    //           desc: '',
    //         },
    //         proxy: {
    //           label: 'OpenAI Proxy',
    //           key: 'proxy',
    //           value: '',
    //           type: TemplateArgType.TextField,
    //           hint: 'Input please',
    //           desc: '',
    //         },
    //         temperature: {
    //           label: 'Temperature',
    //           key: 'temperature',
    //           value: '0.7',
    //           type: TemplateArgType.TextField,
    //           hint: 'Input please',
    //           desc: '',
    //         },
    //       },
    //       '0.0.1',
    //     ),
    //   );
    // }
  }

  async updateLLM(uuid: string, update: UpdateLlmsDto) {
    return await this.dataSource.manager
      .createQueryBuilder()
      .from(LlmsEntity, 'llm')
      .update()
      .set({
        fields: update.fields,
      })
      .where({ uuid })
      .execute();
  }

  async getLLM(uuid: string) {
    console.log('uuid', uuid);
    return await this.dataSource.manager.findOne(LlmsEntity, {
      where: { uuid },
    });
  }
}
