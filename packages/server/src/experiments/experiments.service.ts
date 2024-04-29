import { Injectable, Logger, OnModuleInit } from '@nestjs/common';
import { DataSource } from 'typeorm';
import { ExperimentEntity } from './entities/experiment.entity';
import { ExperimentDto } from './dto/experiment.dto';
import { ToolsService } from '../tools/tools.service';
import { ConversationService } from 'src/chats/conversation.service';
import {
  DeviceType,
  LLMType,
  NodeChild,
  NodeTypes,
  ToolType,
} from '@xpcr/common';
import { ConversationEntity } from '../chats/entities/conversation.entity';
import { User } from '../user/entities/user.entity';
import { MessageEntity, Role } from 'src/chats/entities/message.entity';
import * as path from 'path';
import * as fs from 'fs';

@Injectable()
export class ExperimentsService implements OnModuleInit {
  logger = new Logger(ExperimentsService.name);
  constructor(
    private readonly dataSource: DataSource,
    private readonly toolsService: ToolsService,
    private readonly chatsService: ConversationService,
  ) {}

  async getAllExperiments(user: User) {
    const exps = await this.dataSource.manager.findBy(ExperimentEntity, {
      creator: user.uuid,
    });
    return exps.map((exp) => {
      const expDto = exp as ExperimentDto;
      expDto.creatorName = 'MGI X-Teams';
      return expDto;
    });
  }

  async addExperiment(user: User, isFromInit = false) {
    const nodeChildren: NodeChild[] = [];
    // add LLM (only Azure)
    const allLlmTemps = await this.toolsService.getLLMTemplates();
    let azureLLM = allLlmTemps.find(
      (llmTemp) => llmTemp.llmType === LLMType.Azure,
    );
    if (azureLLM) {
      azureLLM = await this.dataSource.manager.save(azureLLM);
      nodeChildren.push({
        type: azureLLM.type,
        name: 'LLM',
        uuid: azureLLM.uuid,
        functionName: azureLLM.name,
        functionDesc: azureLLM.desc,
        required: true,
      });
    }
    // add OT2
    const allDeviceTemps = await this.toolsService.getDeviceTemplates();
    let ot2 = allDeviceTemps.find(
      (deviceTemp) => deviceTemp.deviceType === DeviceType.OT2,
    );
    if (ot2) {
      ot2 = await this.dataSource.manager.save(ot2);
    }
    // add Google search tool
    const allToolTemps = await this.toolsService.getToolTemplates();
    let googleSearch = allToolTemps.find(
      (toolTemp) => toolTemp.toolType === ToolType.GOOGLE_SEARCH,
    );
    if (googleSearch) {
      googleSearch = await this.dataSource.manager.save(googleSearch);
    }

    // add Wikipedia search tool
    let wikipediaSearch = allToolTemps.find(
      (toolTemp) => toolTemp.toolType === ToolType.WIKIPEDIA_SEARCH,
    );
    if (wikipediaSearch) {
      wikipediaSearch = await this.dataSource.manager.save(wikipediaSearch);
    }

    // add Agent
    let allAgentTemps = await this.toolsService.getAgentTemplates();
    allAgentTemps.forEach((agentTemp) => {
      // add LLM and OT2 to Agent
      agentTemp.children.forEach((child) => {
        if (child.type === NodeTypes.LLMs) {
          child.uuid = azureLLM?.uuid || '';
          child.functionName = azureLLM?.name || '';
          child.functionDesc = azureLLM?.desc || '';
        } else if (child.type === NodeTypes.Devices) {
          child.uuid = ot2?.uuid || '';
          child.functionName = ot2?.name || '';
          child.functionDesc = ot2?.desc || '';
        } else if (child.type === NodeTypes.Tools) {
          if (child.name === 'Google Search') {
            child.uuid = googleSearch?.uuid || '';
            child.functionName = googleSearch?.name || '';
            child.functionDesc = googleSearch?.desc || '';
          } else if (child.name === 'Wikipedia Search') {
            child.uuid = wikipediaSearch?.uuid || '';
            child.functionName = wikipediaSearch?.name || '';
            child.functionDesc = wikipediaSearch?.desc || '';
          }
        }
      });
    });
    allAgentTemps = await this.dataSource.manager.save(allAgentTemps);

    allAgentTemps.forEach((agentTemp) => {
      nodeChildren.push({
        type: agentTemp.type,
        name: agentTemp.name,
        uuid: agentTemp.uuid,
        functionName: agentTemp.name,
        functionDesc: agentTemp.desc,
        required: true,
      });
    });
    const conversationName = isFromInit
      ? `Multiplex PCR(${user.userName})`
      : `New Multiplex PCR`;
    const conversation = await this.dataSource.manager.save(
      new ConversationEntity(
        conversationName,
        'Modify the board layout of OT-2.',
        user.uuid,
      ),
    );
    const experimentName = isFromInit
      ? `Multiplex PCR With OT-2 (${user.userName})`
      : 'New Multiplex PCR With OT-2';
    await this.dataSource.manager.save(
      new ExperimentEntity(
        experimentName,
        nodeChildren,
        'Simultaneously amplifies multiple target DNA sequences using the OT-2 platform for enhanced accuracy and efficiency.',
        user.uuid,
        conversation.uuid,
      ),
    );

    // add intro
    const introText = fs
      .readFileSync(
        path.join(process.cwd(), 'assets', 'plans', 'introduce.json'),
      )
      .toString('utf-8');
    await this.chatsService.createMessage(
      '',
      JSON.parse(introText)['EN'],
      conversation.uuid,
      Role.Assistant,
    );
  }

  async onModuleInit() {
    const users = await this.dataSource.manager.find(User);
    if (users && users.length > 0) {
      for (const user of users) {
        const exists = await this.dataSource.manager.exists(ExperimentEntity, {
          where: { creator: user.uuid },
        });
        if (exists) {
          this.logger.warn(`experiments for ${user.userName} exists`);
          continue;
        }
        await this.initMultiPCR(user);
        await this.initDNALib(user);
      }
    }
  }

  private async initDNALib(user: User) {
    const conversation = await this.dataSource.manager.save(
      new ConversationEntity(
        'DNA Preparation',
        'Report of DNA Preparation.',
        user.uuid,
      ),
    );
    await this.dataSource.manager.save(
      new ExperimentEntity(
        'DNA Library Preparation',
        [],
        'DNA library preparation is a crucial process in genomics research, where fragmented DNA samples are prepared for high-throughp...',
        user.uuid,
        conversation.uuid,
      ),
    );
  }

  private async initMultiPCR(user: User) {
    this.addExperiment(user, true);
  }
}
