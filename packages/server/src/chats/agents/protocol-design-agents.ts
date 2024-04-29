import { BaseAgents } from './base-agents';
import { ChatMessage } from '@azure/openai/types/src';
import { AgentFunctions } from '@xpcr/common';
import { Logger } from '@nestjs/common';
import { Agents } from '../../tools/agents/entities/agents.entity';
import { DataSource } from 'typeorm';
import * as fs from 'fs';
import * as path from 'path';
import { Role } from '../entities/message.entity';

const mockResult = `
  根据你的输入"{{INPUT}}",
  设计出的Protocol如下:
  
\`\`\`python
from opentrons.types import Point
metadata = {
    'protocolName': 'Fault detection demo',
    'author': 'MGI-X',
    'apiLevel': '2.11'
}

def run(ctx):
    # Modules
    mag_mod = ctx.load_module('magnetic module gen2', '1')
    mag_rack = mag_mod.load_labware('biorad_96_wellplate_200ul_pcr') 
    tmp_mod = ctx.load_module('temperature module gen2', '3')
    tmp_rack=tmp_mod.load_labware('biorad_96_wellplate_200ul_pcr') 
    tc_mod = ctx.load_module(module_name='thermocyclerModuleV1')
    tc_rack = tc_mod.load_labware(name='biorad_96_wellplate_200ul_pcr')
    # pipette and tiprack
    tiprack = [ctx.load_labware('opentrons_96_filtertiprack_20ul', slot) for slot in ['2','5']]
    pip20_multi = ctx.load_instrument('p20_multi_gen2', 'left', tip_racks=[*tiprack])
    pip20_single = ctx.load_instrument('p20_single_gen2', 'right', tip_racks=[*tiprack])
    # well plates
    plate_4 = ctx.load_labware('nest_96_wellplate_2ml_deep', '4', 'reg4')
    plate_6 = ctx.load_labware('biorad_96_wellplate_200ul_pcr', '6', 'reg6')
    plate_9 = ctx.load_labware('biorad_96_wellplate_200ul_pcr', '9', 'reg9')
    # protocol
    for well_name in ['A1','B1','C1']:
        pip20_single.transfer(10,
                            mag_rack.wells_by_name()[well_name],
                            tmp_rack.wells_by_name()[well_name],
                            new_tip='always',
                            blow_out=True,
                            blowout_location='destination well')
    for well_name in ['A1','A2','A3']:
        pip20_multi.transfer(10,
                            mag_rack.wells_by_name()[well_name],
                            tmp_rack.wells_by_name()[well_name],
                            new_tip='always',
                            blow_out=True,
                            blowout_location='destination well')
\`\`\`
`;

export class ProtocolDesignAgents extends BaseAgents {
  initMessages: ChatMessage[] = [];
  availableFunctions = {
    [AgentFunctions.DESIGN_PRIMER]: this.designProtocol.bind(this),
  };

  constructor(readonly agent: Agents, readonly dataSource: DataSource) {
    super(agent, dataSource);
    this.needSummarize = false;
  }

  async init(conversationUUID) {
    const protocolDesignAgent = fs.readFileSync(
      path.join(
        process.cwd(),
        'assets',
        'prompts',
        'protocol-design-agent.txt',
      ),
    );
    this.initMessages = [
      {
        role: Role.System,
        content: String(protocolDesignAgent),
      },
    ];
    await super.init(conversationUUID);
  }

  // 引物设计不需要子Agent内部的LLM进行总结，提升反馈的效率
  async *send(
    userInput: string,
    description: string,
  ): AsyncGenerator<any, void, any> {
    // yield* this.mockGenerating(`[${this.agent.name}]\n`);
    let content = '';
    for (const key in this.availableFunctions) {
      const queryFunc = this.availableFunctions[key];
      for await (const msg of queryFunc({
        query: `${userInput} ${description}`,
      })) {
        // 将用户完整的输入传给Seq Agent
        if (msg.role === Role.Assistant) {
          content += msg.content;
        } else {
          yield msg;
        }
      }
    }
    yield {
      role: Role.Assistant,
      content: content,
    };
    await this.saveMessage(`${userInput} ${description}`, Role.User);
    await this.saveMessage(content, Role.Assistant);
  }

  private async *designProtocol({ query }: { query: string }) {
    this.logger.debug(`designProtocol ${query}`);
    let success = true;
    let protocolResult = '';
    let python: string;

    try {
      protocolResult = mockResult
        .replace('{{INPUT}}', query)
        .replace(/^\s+/gm, '');
      this.logger.debug('protocolResult\n' + protocolResult);

      // 提取python代码块
      // const codeBlockRegex = /```(?:\w+)?\s*([\s\S]+?)\s*```/g;
      // const codeBlocks = mockResult.match(codeBlockRegex);
      // if (codeBlocks) {
      //   codeBlocks.forEach((codeBlock, index) => {
      //     this.logger.debug(codeBlock);
      //     python = JSON.parse(codeBlock.replace(/```python|```/g, ''));
      //   });
      // }
    } catch (e) {
      this.logger.error(`designProtocol error: ${e}`, e.stack);
      protocolResult = e.message;
      success = false;
    }

    yield {
      content: protocolResult,
      role: Role.Assistant,
      protocol: python,
    };
  }

  getLogger(): Logger {
    return new Logger(ProtocolDesignAgents.name);
  }
  getInitMessages(): ChatMessage[] {
    return this.initMessages;
  }
  getAvailableFunctions(): {
    [p: string]: (params: object) => AsyncGenerator<any, void, any>;
  } {
    return this.availableFunctions;
  }
}
