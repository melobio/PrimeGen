import { Injectable, Logger, OnModuleInit } from '@nestjs/common';
import { DataSource } from 'typeorm';

@Injectable()
export class SeedAgentsService implements OnModuleInit {
  logger = new Logger(SeedAgentsService.name);
  constructor(private readonly dataSource: DataSource) {}

  // async create() {
  //   const exists = await this.dataSource.manager.exists(Agents);
  //   if (exists) {
  //     this.logger.warn(`Agents exists.`);
  //   } else {
  //     // Seed Agents
  //     await this.dataSource.manager.save(
  //       new Agents(
  //         AGENT_FAULT_NAME,
  //         'Check whether the OT2 or OT3 machine has encountered a runtime error based on the provided parameters.',
  //         'http://127.0.0.1:3000',
  //         [
  //           {
  //             name: CHECK_OT2_STATE,
  //             description:
  //               'Check whether the OT2 machine has encountered a runtime error based on the provided parameters.',
  //             parameters: {
  //               type: 'object',
  //               properties: {
  //                 checkType: {
  //                   type: 'string',
  //                   description: "The type of the check, 'tips' or 'liquid'.",
  //                 },
  //               },
  //             },
  //             required: ['type'],
  //           },
  //         ],
  //         [],
  //         '0.0.1',
  //       ),
  //     );
  //     await this.dataSource.manager.save(
  //       new Agents(
  //         'Google Search Agent',
  //         'Google Search Agent is a tool designed to automate and optimize search queries on the Google platform, enh...',
  //         'http://127.0.0.1:3000',
  //         [],
  //         [],
  //         '0.0.1',
  //       ),
  //     );
  //     await this.dataSource.manager.save(
  //       new Agents(
  //         'Documents Search Agent',
  //         'Documents Search Agent is a specialized tool for efficiently scanning and retrieving specific content from....',
  //         'http://127.0.0.1:3000',
  //         [],
  //         [],
  //         '0.0.1',
  //       ),
  //     );
  //     await this.dataSource.manager.save(
  //       new Agents(
  //         'Protocol Optimization Agent',
  //         'Protocol Optimization Agent is a sophisticated tool that dynamically tunes communication protocols to e....',
  //         'http://127.0.0.1:3000',
  //         [],
  //         [],
  //         '0.0.1',
  //       ),
  //     );
  //     await this.dataSource.manager.save(
  //       new Agents(
  //         AGENT_CODE_EXECUTION_NAME,
  //         'Code Execution Agent is a software component designed to execute the protocol on OT2 and return the results.',
  //         'http://127.0.0.1:31950',
  //         [
  //           {
  //             name: EXECUTE_OT2_PROTOCOL,
  //             description:
  //               'Execute the protocol on OT2 and return the results.',
  //             parameters: {
  //               type: 'object',
  //               properties: {
  //                 protocol: {
  //                   type: 'string',
  //                   description: 'The protocol to be executed.',
  //                 },
  //               },
  //             },
  //             required: ['protocol'],
  //           },
  //         ],
  //         [],
  //         '0.0.1',
  //       ),
  //     );
  //     await this.dataSource.manager.save(
  //       new Agents(
  //         'Experiment Designer Agent',
  //         'Experiment Designer Agent is a digital assistant crafted to aid researchers in planning, structuring, and executing s....',
  //         'http://127.0.0.1:3000',
  //         [],
  //         [],
  //         '0.0.1',
  //       ),
  //     );
  //     await this.dataSource.manager.save(
  //       new Agents(
  //         'Active Learning Agent',
  //         'An intelligent system designed to adaptively select the most informative data samples, enhancing the learning....',
  //         'http://127.0.0.1:3000',
  //         [],
  //         [],
  //         '0.0.1',
  //       ),
  //     );
  //   }
  // }

  async onModuleInit() {
    // await this.create();
  }
}
