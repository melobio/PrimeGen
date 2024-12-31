import { Logger } from '@nestjs/common';
import * as fs from 'fs';
import * as path from 'path';
import {
  AzureKeyCredential,
  ChatCompletions,
  ChatMessage,
  FunctionCall,
  FunctionDefinition,
  OpenAIClient,
} from '@azure/openai';
import { OpenAI } from 'openai';
import { DataSource, IntegerType } from 'typeorm';
import { ExperimentEntity } from '../../experiments/entities/experiment.entity';
import {
  AgentFunctionByType,
  AgentFunctions,
  AgentType,
  LLMFields,
  LLMType,
  NodeTypes,
} from '@xpcr/common';
import {
  CancerOptionInfo,
  GeneticDisorderInfo,
  NCBIFunctionType,
  NCBIResultContent,
  PathogenDrugInfo,
  ProteinMutationInfo,
  SnpPrimerDesignInfo,
  SpeciesIdentificationInfo,
} from '../agents/ext/ncbi-search';
import { ProtocolDesignInfo } from '../agents/ext/protocol-design';
import { LlmsEntity } from '../../tools/llms/entities/llms.entity';
import { Agents } from '../../tools/agents/entities/agents.entity';
import { FaultAgents } from './fault-agents';
import { CodeExecutionAgent } from './code-execution-agent';
import { BaseAgents } from './base-agents';
import { Socket } from 'socket.io';
import {
  MessageEntity,
  MiddleMessageType,
  Role,
} from '../entities/message.entity';
import * as _ from 'lodash';
import {
  TEXT_ANSWER_CREATE,
  TEXT_ANSWER_DONE,
  TEXT_ANSWER_GENERATING,
  UPLOAD_PRIMER_EXCEL,
  INITIATIVE_START_PRIMER_DESIGN,
  RESTART_CONVERSATION,
  UPDATE_CONVERSATIONS,
  UPDATE_MESSAGE,
  SUMMARY_MESSAGE,
  INITIATIVE_START_CODE_EXECUTION,
} from '../conversation.constant';
import { InternetSearchAgents } from './internet-search-agents';
import { SequenceSearchAgents } from './sequence-search-agents';
import { ConversationService } from '../conversation.service';
import { PlanStep, Planner } from './ext/planner';
import { ConversationEntity } from '../entities/conversation.entity';
import { PrimerDesignAgents } from './primer-design-agents';
import { ProtocolDesignAgents } from './protocol-design-agents';
import { CacheService } from 'src/cache/cache.service';
import { PrimerDesign, PrimerResultContent } from './ext/primer-design';
import { AgentMessageEntity } from '../entities/agent-message.entity';
interface jasonType {
  output: string;
  achieved: boolean;
  thought: string;
  objective?: string;
  wantTips?: string;
}
function removeLeadingSpaces(str) {
  // use regex to remove leading spaces
  return str.replace(/^\s+/gm, '');
}

function removeLineBreaks(str) {
  // use regex to remove line breaks
  return str.replace(/\n/g, '');
}

function checkGeneFile(inputString: string) {
  // use regex to match file name
  const fileRegex = /\.(fasta|fastq|fna|csv)\b/;
  return fileRegex.test(inputString);
}
// use regex to match gene info if purpose is not achieved
function checkStringForFileOrGeneInfo(inputString: string) {
  // regex to match gene info file
  const geneInfoRegex1 =
    /\b(?=.*chrom)(?=.*start)(?=.*end)(?=.*gene)(?=.*type)(?=.*HGNC)\b/;
  const geneInfoRegex2 =
    /\b(?=.*Gene name)(?=.*Chromosome)(?=.*Gene start)(?=.*Gene end)(?=.*Start position in gene)(?=.*End position in gene)\b/;
  // check if input string matches file name or gene info fields
  if (
    checkGeneFile(inputString) ||
    geneInfoRegex1.test(inputString) ||
    geneInfoRegex2.test(inputString)
  ) {
    return true;
  } else {
    return false;
  }
}

export class RouterAgents {
  logger = new Logger(RouterAgents.name);
  planner: Planner;
  openAIClient: OpenAIClient;
  GPTClient: OpenAI;
  llm_type: LLMType;
  initMessages: ChatMessage[] = [];
  // FUNCTION to sub agent
  functionToAgents: { [key: string]: BaseAgents } = {};
  examples: { [step: string]: string } = {};
  checkPurpose: string;
  next: string;
  summary: string;
  cvs: ConversationEntity;
  constructor(
    private readonly dataSource: DataSource,
    private readonly cacheService: CacheService,
    private readonly experiment: ExperimentEntity,
    private readonly chatsService: ConversationService,
    private readonly conversationUUID: string,
  ) {}

  async updateTokens(usage) {
    if (!usage) {
      return;
    }
    this.cvs.completionTokens += usage.completionTokens;
    this.cvs.promptTokens += usage.promptTokens;
    this.cvs.totalTokens += usage.totalTokens;
    console.log(
      `[${this.conversationUUID}] completionTokens: ${this.cvs.completionTokens}`,
    );
    console.log(
      `[${this.conversationUUID}] promptTokens: ${this.cvs.promptTokens}`,
    );
    console.log(
      `[${this.conversationUUID}] totalTokens: ${this.cvs.totalTokens}`,
    );
    this.dataSource.manager.update(ConversationEntity, this.cvs.id, this.cvs);
  }

  // update step result
  updateStepResult(stepName: string, result: string) {
    this.cvs.stepsResult[stepName] = result;
    this.dataSource.manager
      .update(ConversationEntity, this.cvs.id, this.cvs)
      .then();
  }

  async init() {
    this.cvs = await this.dataSource.manager.findOneBy<ConversationEntity>(
      ConversationEntity,
      {
        uuid: this.conversationUUID,
      },
    );
    if (!this.cvs.stepsResult) {
      this.cvs.stepsResult = {};
    }

    // read check-purpose template content
    this.checkPurpose = fs
      .readFileSync(
        path.join(process.cwd(), 'assets', 'plans', 'check-purpose.txt'),
      )
      .toString('utf-8');
    // this.logger.debug(`checkPurpose: ${String(this.checkPurpose)}`);

    // read next template content
    this.next = fs
      .readFileSync(path.join(process.cwd(), 'assets', 'plans', 'next.txt'))
      .toString('utf-8');
    // this.logger.debug(`next: ${String(this.next)}`);

    // read summary template content
    this.summary = fs
      .readFileSync(path.join(process.cwd(), 'assets', 'plans', 'summary.txt'))
      .toString('utf-8');
    // this.logger.debug(`summary: ${String(this.summary)}`);

    // read PCR experiment preset data, initialize planner, and set cvs's currentStep
    this.planner = Planner.parse(
      fs
        .readFileSync(path.join(process.cwd(), 'assets', 'plans', 'plan.json'))
        .toString('utf-8'),
      this.cvs,
      this.dataSource,
    );
    // this.planner.setCurrentStep(cvs.currentStep);
    const routerAgent = fs.readFileSync(
      path.join(process.cwd(), 'assets', 'prompts', 'router-agent.txt'),
    );
    // this.logger.debug(`routerAgent: ${String(routerAgent)}`);
    this.initMessages = [
      {
        role: Role.System,
        content: String(routerAgent).replace(
          '{{STEPS}}',
          JSON.stringify(this.planner.steps),
        ),
      },
    ];
    // this.logger.debug('initMessages:', this.initMessages);
    await this._initOpenAI();
    await this._loadAgents();
  }

  private async _initOpenAI() {
    let apiKey = '';
    let apiBase = '';
    let apiEngine = '';
    const llmChildren = this.experiment.nodeChildren.filter(
      (nodeChild) => nodeChild.type === NodeTypes.LLMs,
    );
    // init azure llm
    for (const llmChild of llmChildren) {
      const llm = await this.dataSource.manager.findOne(LlmsEntity, {
        where: { uuid: llmChild.uuid },
      });
      if (llm && llm.llmType === LLMType.Azure) {
        llm.fields?.forEach((field) => {
          if (field.key === LLMFields.API_KEY) {
            apiKey = field.value;
          } else if (field.key === LLMFields.API_BASE) {
            apiBase = field.value;
          } else if (field.key === LLMFields.API_ENGINE) {
            apiEngine = field.value;
          }
        });
      }
    }
    if (!apiKey || !apiBase || !apiEngine) {
      apiKey = process.env.OPENAI_API_KEY;
      apiBase = process.env.OPENAI_API_BASE;
      apiEngine = process.env.OPENAI_API_ENGINE;
      this.logger.warn(`LLM config not found, use default config.`);
    }
    this.logger.log(`open ai key: ${apiKey}`);
    this.logger.log(`open ai base: ${apiBase}`);
    this.logger.log(`open ai engine: ${apiEngine}`);
    if (process.env.LLM_TYPE.toLowerCase() == 'azure') {
      this.llm_type = LLMType.Azure;
      this.openAIClient = new OpenAIClient(
        apiBase,
        new AzureKeyCredential(apiKey),
      );
    } else if (process.env.LLM_TYPE.toLowerCase() == 'openai') {
      this.llm_type = LLMType.OpenAI;
      this.GPTClient = new OpenAI({
        baseURL: apiBase,
        apiKey: apiKey,
      });
    }
    this.logger.log(`init openai success.`);
  }

  private async _loadAgents() {
    const agentChildren = this.experiment.nodeChildren.filter(
      (nodeChild) => nodeChild.type === NodeTypes.Agents,
    );
    for (const agentChild of agentChildren) {
      const agent = await this.dataSource.manager.findOne(Agents, {
        where: { uuid: agentChild.uuid },
      });
      if (agent) {
        const funcName = AgentFunctionByType[agent.agentType];
        let newAgent: BaseAgents;
        switch (agent.agentType) {
          case AgentType.FAULT:
            newAgent = new FaultAgents(
              agent,
              this.dataSource,
              this.chatsService,
            );
            break;
          case AgentType.CODE_EXECUTION:
            newAgent = new CodeExecutionAgent(
              agent,
              this.dataSource,
              this.cacheService,
              this.chatsService,
            );
            break;
          case AgentType.INTERNET_SEARCH:
            newAgent = new InternetSearchAgents(
              agent,
              this.dataSource,
              this.cacheService,
            );
            break;
          case AgentType.SEQUENCE_SEARCH:
            newAgent = new SequenceSearchAgents(
              agent,
              this.dataSource,
              this.cacheService,
              this.chatsService,
            );
            break;
          case AgentType.PRIMER_DESIGN:
            newAgent = new PrimerDesignAgents(
              agent,
              this.dataSource,
              this.cacheService,
              this.chatsService,
            );
            break;
          case AgentType.PROTOCOL_DESIGN:
            newAgent = new ProtocolDesignAgents(
              agent,
              this.dataSource,
              this.cacheService,
              this.chatsService,
            );
            break;
        }
        if (newAgent) {
          await newAgent.init(this.conversationUUID);
          this.functionToAgents[funcName] = newAgent;
        }
        // console.log('AgentFunctionByType', AgentFunctionByType);
        this.logger.log(`load agent [${funcName}] -> [${agent.agentType}]`);
      }
    }
  }

  private _getFunctions(): FunctionDefinition[] {
    const current = this.planner.getCurrentStep();
    this.logger.debug('current use agent', current?.useAgents);
    const funcDefs: FunctionDefinition[] = [];
    for (const funcName in this.functionToAgents) {
      const agent: BaseAgents = this.functionToAgents[funcName];

      if (
        !current?.useAgents ||
        !current?.useAgents.includes(agent.agent.agentType)
      ) {
        continue;
      }
      funcDefs.push({
        name: funcName,
        description: agent.agent.desc,
        parameters: {
          type: 'object',
          properties: {
            description: {
              type: 'string',
              description: `The original input of the user, don't summarize and change`,
            },
          },
        },
      });
    }
    // console.log(funcDefs);
    return funcDefs;
  }

  // finish the first step of the STEP, need to prompt the user the next step
  async showSummary(client: Socket, history: MessageEntity[]) {
    const nextStep = this.planner.getNextStep();

    const assistantMessage = await this.chatsService.createMessage(
      SUMMARY_MESSAGE,
      '',
      this.conversationUUID,
      Role.Assistant,
    );
    // send to client
    this.chatsService.emitToCliet(client, TEXT_ANSWER_CREATE, {
      success: true,
      data: assistantMessage,
      finishReason: null,
      conversationUUID: this.conversationUUID,
    });
    let content = `
      #known information: ${JSON.stringify(this.cvs.stepsResult)}
      #Objective: Quickly assist in planning the next steps of the experiment for the user.
      #Reflection: Integrate the completed experimental process with the upcoming experimental process to help the user plan the next steps of the experiment.
      #Output: first answer The next step of the experimental processStarting from some connective words, such as "according to your requirements", "according to the experimental requirements" ...;
      second The final output is based on the name of the next step (${
        nextStep.stepName
      }), Combined with the context; Introduce you to the user to assign ${
      nextStep.useAgents
    }agent will help user to achieve the next goal;
    third introduce begin ${nextStep.stepName},
    Please use English to reply
      #Note: The output language is consistent with the language entered by the user.
      #Rule
      * Please keep your responses as concise as possible
      * Please use English to reply.
    `;
    content = removeLeadingSpaces(content);
    this.logger.debug(`content: ${content}`);

    const currentStepIndex = this.planner.steps.indexOf(
      this.planner.getCurrentStep(),
    );
    const doneSteps = this.planner.steps.slice(0, currentStepIndex + 1);
    const leftSteps = this.planner.steps.slice(currentStepIndex + 1);

    const system = this.summary
      .replace('{{ALL_STEPS}}', JSON.stringify(this.planner.steps))
      .replace('{{DONE_STEPS}}', JSON.stringify(doneSteps))
      .replace('{{NEXT_STEPS}}', JSON.stringify(leftSteps))
      .replace(
        '{{PURPOSE}}',
        this.cvs.stepsResult[this.planner.steps[0].stepName],
      );
    const startTime = performance.now();
    let events = null;
    if (this.llm_type === LLMType.Azure) {
      const messages = [
        { role: Role.System, content: system },
        { role: Role.User, content },
      ];
      this.logger.debug(`summary messages:${JSON.stringify(messages)}`);
      events = await this.openAIClient.listChatCompletions(
        process.env.OPENAI_API_ENGINE,
        messages,
      );
      this.logger.debug(`this.openAIClient.events: ${JSON.stringify(events)}`);
      this.updateTokens(events.usage);
      await this.handleStreamMessage(client, assistantMessage, events);
    } else if (this.llm_type === LLMType.OpenAI) {
      const messages: OpenAI.ChatCompletionMessageParam[] = [
        { role: Role.System, content: system, name: 'system' },
        { role: Role.User, content, name: 'user' },
      ];
      const params: OpenAI.ChatCompletionCreateParams = {
        model: process.env.OPENAI_API_ENGINE,
        messages: messages,
        stream: true,
      };
      events = await this.GPTClient.chat.completions.create(params);
      await this.handleStreamMessage(client, assistantMessage, events);
    }

    const endTime = performance.now();
    const executionTime = ((endTime - startTime) / 1000).toFixed(2);
    this.logger.warn(`Execution summary time: ${executionTime} seconds`);
  }

  async checkAchieved(client: Socket, input?: string) {
    this.logger.debug(`check whether achieved`);
    const currentStep = this.planner.getCurrentStep();
    const isConfirmExperObj =
      currentStep.stepName == this.planner.steps[0].stepName;
    const isDeterminingDNA =
      currentStep.stepName == this.planner.steps[1].stepName;

    const history = await this.chatsService.findAllMessages(
      this.conversationUUID,
    );
    // remove the opening remarks
    history.splice(0, 1);
    let json: jasonType = {
      output: '',
      achieved: false,
      thought: '',
      objective: '',
      wantTips: '',
    };
    let cvs = history.map((item) => `${item.role}: ${item.content}`).join('\n');
    if (isConfirmExperObj) {
      cvs = history
        .filter((item) => !item.agentType)
        .map((item) => `${item.role}: ${item.content}`)
        .join('\n');
    }
    this.logger.warn(`checkAchieved cvs:${JSON.stringify(cvs)}`);
    // extract the experimental purpose
    let content = `
    #Purpose:  Careful and comprehensive confirm: ${
      currentStep.purpose
    } from the dialogue;\n
    ${
      currentStep.purposeExample
        ? `Examples:\n ${currentStep.purposeExample}`
        : ''
    }
    #Think: Analyze the dialogue to determine whether the assistant can get the purpose from the dialogue. Please don't guess or suggest.\n
    #Output：Please quickly output a JSON object with the following keys (do not start or end with \`\`\`json):\n
    - output //  String type, the purpose summarized from the dialogue. It needs to include the actual results.\n
    - achieved // Boolean type, whether the purpose was achieved.\n
    - thought // String type, the thinking process.\n
    - wantTips // String type, when the purpose is not achieved, provide a prompt to help achieve the purpose. Explain: "So I can better help you complete ${
      currentStep.purpose
    });\n
    ${
      isConfirmExperObj
        ? `- objective // String type, the experimental purpose (composed of experiment type and experimental subject).\n`
        : ''
    }
    # Rule:
    * Please use English to reply.`;
    content = removeLeadingSpaces(content);
    this.logger.debug(`content: ${content}`);
    let system = '';
    if (currentStep.checkWithContext) {
      system = this.checkPurpose.replace('{{CVS}}', cvs);
    } else {
      const agentCvs = history
        .filter(
          (item) =>
            item.role == Role.Assistant &&
            currentStep.useAgents.includes(item.agentType),
        )
        .map((item) => {
          let optStr = '';
          if (item?.optionInfo && item.optionInfo.state == 'stop') {
            optStr = JSON.stringify({ ...item.optionInfo, operations: [] });
          }
          return `${item.role}: ${removeLineBreaks(item.content)};${optStr}`;
        })
        .join('\n');
      system = this.checkPurpose.replace('{{CVS}}', agentCvs);
    }
    this.logger.debug(`system: ${system}`);

    const startTime = performance.now();
    let events = null;
    if (this.llm_type == LLMType.Azure) {
      const messages = [
        { role: Role.System, content: system },
        { role: Role.User, content },
      ];
      this.logger.warn(`checkAchieved messages:${JSON.stringify(messages)}`);
      events = await this.openAIClient.getChatCompletions(
        process.env.OPENAI_API_ENGINE,
        messages,
      );
      this.logger.log(`this.openAIClient.events: ${JSON.stringify(events)}`);
      this.updateTokens(events.usage);
    } else {
      const params: OpenAI.ChatCompletionCreateParams = {
        model: process.env.OPENAI_API_ENGINE,
        messages: [
          { role: 'system', content: system, name: 'system' },
          { role: 'user', content: content, name: 'user' },
        ],
        stream: false,
      };
      this.logger.warn(`checkAchieved messages: ${JSON.stringify(params)}`);
      events = await this.GPTClient.chat.completions.create(params);
    }
    const endTime = performance.now();
    const executionTime = ((endTime - startTime) / 1000).toFixed(2);
    this.logger.warn(`Execution checkAchieved time: ${executionTime} seconds`);

    if (events.choices && events.choices.length > 0) {
      const choice = events.choices[0];

      if (choice.message) {
        this.logger.warn(`check message: ${choice.message.content}`);

        try {
          json = JSON.parse(choice.message.content);
        } catch (e) {
          this.logger.warn(`JSON parse error: ${e.message}`);
          this.logger.warn(`Use reg match instead.`);
          const codeBlockRegex = /```(?:\w+)?\s*([\s\S]+?)\s*```/g;
          const codeBlocks = choice.message.content.match(codeBlockRegex);
          if (codeBlocks) {
            codeBlocks.forEach((codeBlock, index) => {
              this.logger.debug(codeBlock);
              json = JSON.parse(codeBlock.replace(/```json|```/g, ''));
            });
          }
        }

        this.logger.debug(`json: ${JSON.stringify(json)}`);
      }
    }
    if (json.achieved) {
      this.updateStepResult(currentStep.stepName, json.output);
      // await this.showSummary(client, history);
      this.planner.goNextStep();
      if (isConfirmExperObj && json.objective) {
        // update conversation name
        this.cvs.name = json.objective;
        this.dataSource.manager.update(
          ConversationEntity,
          this.cvs.id,
          this.cvs,
        );
        // update experiment name
        const experiment = await this.dataSource.manager.findOne(
          ExperimentEntity,
          {
            where: {
              conversationUUID: this.cvs.uuid,
            },
          },
        );
        experiment.name = `${json.objective} With OT-2`;
        this.dataSource.manager.update(
          ExperimentEntity,
          experiment.id,
          experiment,
        );
        client.emit(UPDATE_CONVERSATIONS);
        // after knowing the purpose, send to search agent
        await this.sendToSearch(client, input);
      }
    }
    return json.achieved;
  }

  async handleNcbiSearchStop(
    client: Socket,
    responses: NCBIResultContent['responses'],
  ) {
    if (responses.stage > 2) {
      const responsesStr = JSON.stringify(responses);
      if (checkStringForFileOrGeneInfo(responsesStr)) {
        const history = await this.chatsService.findAllMessages(
          this.conversationUUID,
        );
        const currentStep = this.planner.getCurrentStep();
        this.updateStepResult(currentStep.stepName, responsesStr);
        // 24.12.19 remove Reflection
        // await this.showSummary(client, history);
        this.planner.goNextStep();
      } else {
        const achieved = await this.checkAchieved(client, responsesStr);
        if (!achieved) {
          this.logger.debug('=====ncbi-search=checkAchieved=fail==');
        }
      }
      // next step: primer-design, create primer-design tips
      await this.createPrimerDesignTips(client);
    } else {
      this.logger.debug('sendTextToSearch==>search end, start normal chat');
    }
  }

  async handlePrimerDesignStop(
    client: Socket,
    msg: {
      optionInfo: SnpPrimerDesignInfo;
      content: string;
      role: Role;
    },
  ) {
    if (msg.optionInfo.stage > 2) {
      if (checkGeneFile(msg.content)) {
        const history = await this.chatsService.findAllMessages(
          this.conversationUUID,
        );
        const currentStep = this.planner.getCurrentStep();
        this.updateStepResult(currentStep.stepName, msg.content);
        await this.showSummary(client, history);
        this.planner.goNextStep();
      } else {
        const achieved = await this.checkAchieved(client, msg.content);
        if (!achieved) {
          this.logger.debug('=====primer-design=checkAchieved=fail==');
        }
      }
    } else {
      this.logger.debug(
        'sendTextToSearch ==>primer-design end, start normal chat',
      );
    }
  }
  // create start PrimerDesign tip message
  async createPrimerDesignTips(client: Socket) {
    await this.chatsService.createMessage(
      INITIATIVE_START_PRIMER_DESIGN,
      "We have now completed the 'Determining DNA Sequence' step. Would you like to start the primer design process now?",
      this.conversationUUID,
      Role.Assistant,
    );
    client.emit(UPDATE_MESSAGE, {
      conversationUUID: this.conversationUUID,
    });
  }

  async sendToCodeExecution(client: Socket, msg: string, optionInfo?: any) {
    const previewFunctionCall = {
      name: AgentFunctionByType[AgentType.CODE_EXECUTION],
      arguments: `{"description":"${msg}"}`,
    };
    const agent: BaseAgents = this.functionToAgents[previewFunctionCall.name];
    const agentMessage = await this.chatsService.createMessage(
      '',
      '',
      this.conversationUUID,
      Role.Assistant,
      agent.agent.agentType,
    );
    this.chatsService.emitToCliet(client, TEXT_ANSWER_CREATE, {
      success: true,
      data: agentMessage,
      finishReason: null,
      conversationUUID: this.conversationUUID,
    });
    for await (const response of this._callFunction(
      previewFunctionCall,
      msg,
      optionInfo,
    )) {
      if (response.role == Role.Assistant) {
        // this.logger.log(
        //   `CodeExecution response-optionInfo: ${JSON.stringify(
        //     response.optionInfo,
        //   )}`,
        // );
        if (response?.optionInfo) {
          agentMessage.optionInfo = response.optionInfo;
        }
        // string content
        agentMessage.content = response.content;
        this.chatsService.emitToCliet(client, TEXT_ANSWER_GENERATING, {
          success: true,
          data: agentMessage,
          finishReason: null,
          conversationUUID: this.conversationUUID,
        });
        if (!response?.optionInfo?.data?.generating) {
          // need user to choose what to do next
          this.chatsService.emitToCliet(client, TEXT_ANSWER_GENERATING, {
            success: true,
            data: agentMessage,
            finishReason: 'stop',
            conversationUUID: this.conversationUUID,
          });
          client.emit(UPDATE_MESSAGE, {
            conversationUUID: this.conversationUUID,
          });
          client.emit(TEXT_ANSWER_DONE, {
            conversationUUID: this.conversationUUID,
          });
        }
        // store in database
        await this.dataSource.manager.update(
          MessageEntity,
          agentMessage.id,
          agentMessage,
        );
        // experiment result
        if (response?.optionInfo?.data.state == 'stop') {
          this.chatsService.emitToCliet(client, TEXT_ANSWER_GENERATING, {
            success: true,
            data: agentMessage,
            finishReason: 'stop',
            conversationUUID: this.conversationUUID,
          });
          const responsesStr = JSON.stringify(response);
          const history = await this.chatsService.findAllMessages(
            this.conversationUUID,
          );
          const currentStep = this.planner.getCurrentStep();
          this.updateStepResult(currentStep.stepName, responsesStr);
          const achieved = await this.checkAchieved(client, responsesStr);
          if (!achieved) {
            this.logger.debug('=====code-execution=checkAchieved=fail==');
          } else {
            await this.showSummary(client, history);
            this.planner.goNextStep();
          }
        }
      }
    }
  }

  async sendText(
    client: Socket,
    userMessageContent: string,
    history: MessageEntity[],
  ) {
    let content = removeLeadingSpaces(userMessageContent);
    // replace " with '
    content = content.replace(/"/g, "'");

    const restart = await this.checkIfRestart(content);
    if (restart) {
      // restart the experiment
      client.emit(RESTART_CONVERSATION, {
        conversationUUID: this.conversationUUID,
      });
      this.cvs.totalTokens =
        this.cvs.completionTokens =
        this.cvs.promptTokens =
          0;
      return;
    }

    this.logger.log(`content: ${content}`);

    const chatHistory = history.map<ChatMessage>((item) => {
      return { role: item.role, content: item.content };
    });
    chatHistory.splice(0, 1);

    const currentStep = this.planner.getCurrentStep();
    this.logger.debug('currentStep====>', JSON.stringify(currentStep));
    const currentStepIndex = this.planner.steps.indexOf(currentStep);
    this.logger.debug('currentStepIndex====>', currentStepIndex);
    const doneSteps =
      currentStepIndex !== 0
        ? this.planner.steps.slice(0, currentStepIndex + 1)
        : [];
    this.logger.debug('doneSteps====>', JSON.stringify(doneSteps));
    this.initMessages.forEach((item, index) => {
      if (item.role == Role.System) {
        let embeddedContent = this.initMessages[index].content;
        embeddedContent = embeddedContent.replace(
          '{{DONE_STEPS}}',
          JSON.stringify(doneSteps),
        );
        embeddedContent = embeddedContent.replace(
          '{{CURRENT_STEP}}',
          JSON.stringify(currentStep),
        );
        this.initMessages[index].content = embeddedContent;
      }
    });
    const messages = [...this.initMessages, ..._.slice(chatHistory, -50)];
    this.logger.log(`messages: ${JSON.stringify(messages)}`);
    const isDeterminingDNA =
      currentStep.stepName == this.planner.steps[1].stepName;
    const isPrimerDesign =
      currentStep.stepName == this.planner.steps[2].stepName;
    // stepName == 'Protocol Design'
    const isProtocolDesign =
      currentStep.stepName == this.planner.steps[3].stepName;
    // stepName == 'Code Execution'
    const isCodeExecution =
      currentStep.stepName == this.planner.steps[4].stepName;
    const searchAgentMessages = history.filter(
      (item) =>
        item.role == Role.Assistant &&
        item.agentType == AgentType.SEQUENCE_SEARCH,
    );
    const primerAgentMessages = history.filter(
      (item) =>
        item.role == Role.Assistant &&
        item.agentType == AgentType.PRIMER_DESIGN,
    );
    const protocolDesignMessages = history.filter(
      (item) =>
        item.role == Role.Assistant &&
        item.agentType == AgentType.PROTOCOL_DESIGN,
    );
    const codeExecutionMessages = history.filter(
      (item) => item.agentType == AgentType.CODE_EXECUTION,
    );
    if (isDeterminingDNA && searchAgentMessages.length > 0) {
      // enter `determine gene` step, invoke Search agent, ignore planner
      await this.sendToSearch(client, content);
    } else if (isPrimerDesign && primerAgentMessages.length > 0) {
      // enter `primer design` step, invoke PrimerDesign agent，ignore planner
      await this.sendToPrimerDesign(client, content);
    } else if (isProtocolDesign && protocolDesignMessages.length > 0) {
      // enter `Protocol Design` step, invoke ProtocolDesign agent, ignore planner
      await this.sendToProtocolDesign(client, content);
    } else if (isCodeExecution) {
      // enter `Code Execution` step, invoke CodeExecution agent, ignore planner
      await this.sendToCodeExecution(client, content);
    } else {
      this.logger.debug('=normal_talk=======>');
      const functions = this._getFunctions();
      this.logger.log(`functions: ${JSON.stringify(functions)}`);
      const assistantMessage = await this.chatsService.createMessage(
        '',
        '',
        this.conversationUUID,
        Role.Assistant,
      );
      this.chatsService.emitToCliet(client, TEXT_ANSWER_CREATE, {
        success: true,
        data: assistantMessage,
        finishReason: null,
        conversationUUID: this.conversationUUID,
      });
      // choose LLM according to llm_type
      let events = null;
      if (this.llm_type === LLMType.Azure) {
        events = await this.openAIClient.listChatCompletions(
          process.env.OPENAI_API_ENGINE,
          messages,
          functions.length > 0
            ? {
                functionCall: 'auto',
                functions: functions,
              }
            : null,
        );
        this.logger.log(`this.openAIClient.events: ${JSON.stringify(events)}`);
        this.updateTokens(events.usage);
      } else if (this.llm_type === LLMType.OpenAI) {
        const chatCompletionMessages: OpenAI.ChatCompletionMessageParam[] =
          messages.map((message) => {
            if ('name' in message) {
              return message;
            } else {
              return { ...message, name: '' }; // or some other default value
            }
          }) as OpenAI.ChatCompletionMessageParam[];
        const params: OpenAI.ChatCompletionCreateParams = {
          model: process.env.OPENAI_API_ENGINE,
          messages: chatCompletionMessages,
          function_call: 'auto',
          functions: functions,
          stream: true,
        };
        events = await this.GPTClient.chat.completions.create(params);
      } else {
        events = null;
      }

      // try catch get openai result
      try {
        await this._handleChatCompletions(
          client,
          content,
          assistantMessage,
          messages,
          events,
        );
      } catch (e) {
        assistantMessage.content = e.message;
        // update database
        await this.dataSource.manager.update(
          MessageEntity,
          assistantMessage.id,
          assistantMessage,
        );
        this.chatsService.emitToCliet(client, TEXT_ANSWER_GENERATING, {
          success: true,
          data: assistantMessage,
          finishReason: null,
          conversationUUID: this.conversationUUID,
        });

        this.chatsService.emitToCliet(client, TEXT_ANSWER_GENERATING, {
          success: true,
          data: assistantMessage,
          finishReason: 'stop',
          conversationUUID: this.conversationUUID,
        });
        throw e;
      }
      await this.checkAchieved(client, content);
    }
  }

  sendToSearch = async (
    client: Socket,
    input: string,
    optionInfo?:
      | CancerOptionInfo
      | GeneticDisorderInfo
      | PathogenDrugInfo
      | SpeciesIdentificationInfo
      | ProteinMutationInfo,
  ) => {
    if (this.chatsService.ncbiSearchs.has(this.conversationUUID)) {
      const previewFunctionCall = {
        name: AgentFunctionByType[AgentType.SEQUENCE_SEARCH],
        arguments: `{"description":"${input}"}`,
      };
      if (!optionInfo) {
        const history = await this.chatsService.findAllMessages(
          this.conversationUUID,
        );
        const searchAgentMessage = history.filter(
          (item) =>
            item.role == Role.Assistant &&
            item.agentType == AgentType.SEQUENCE_SEARCH &&
            item.optionInfo,
        );
        const stageArr = searchAgentMessage
          .map((optItem) => optItem.optionInfo.stage)
          .sort((a, b) => a - b);
        const lastOptionInfo =
          searchAgentMessage[searchAgentMessage.length - 1]?.optionInfo;
        const stage = stageArr.length > 0 ? stageArr.pop() : 1;
        const search_type: NCBIFunctionType = lastOptionInfo?.search_type;
        const protein_mutation_dict =
          lastOptionInfo?.protein_mutation_dict || {};
        const origin_data = lastOptionInfo?.data ?? {};
        optionInfo = {
          stage,
          search_type,
          protein_mutation_dict,
          data: origin_data,
        };
        this.logger.log(`sendToSearch construct ${JSON.stringify(optionInfo)}`);
      }
      const chatHistory = await this.chatsService.findAllMessages(
        this.conversationUUID,
      );
      const chatHistoryArr = chatHistory.map<ChatMessage>((item) => {
        return { role: item.role, content: item.content };
      });
      const operationsObj = {};
      if (optionInfo?.operations) {
        optionInfo.operations.forEach((item) => {
          operationsObj[item.key] = item.value;
        });
      }
      const sumitOption = optionInfo
        ? {
            ...optionInfo,
            ...operationsObj,
            conversation: chatHistoryArr,
          }
        : { ...operationsObj, conversation: chatHistoryArr };
      const agent: BaseAgents = this.functionToAgents[previewFunctionCall.name];
      const agentMessage = await this.chatsService.createMessage(
        '',
        '',
        this.conversationUUID,
        Role.Assistant,
        agent.agent.agentType,
      );
      this.chatsService.emitToCliet(client, TEXT_ANSWER_CREATE, {
        success: true,
        data: agentMessage,
        finishReason: null,
        conversationUUID: this.conversationUUID,
      });
      for await (const msg of this._callFunction(
        previewFunctionCall,
        input,
        sumitOption,
      )) {
        if (msg.role === Role.Assistant) {
          if (msg?.optionInfo) {
            agentMessage.optionInfo = msg.optionInfo;
          }

          agentMessage.content += msg.content;
          this.chatsService.emitToCliet(client, TEXT_ANSWER_GENERATING, {
            success: true,
            data: agentMessage,
            finishReason: null,
            conversationUUID: this.conversationUUID,
          });
          await this.dataSource.manager.update(
            MessageEntity,
            agentMessage.id,
            agentMessage,
          );
          this.chatsService.emitToCliet(client, TEXT_ANSWER_GENERATING, {
            success: true,
            data: agentMessage,
            finishReason: 'stop',
            conversationUUID: this.conversationUUID,
          });
          if (msg?.optionInfo?.state == 'stop') {
            await this.handleNcbiSearchStop(client, msg.optionInfo);
          }
        }
      }
    }
  };

  sendToPrimerDesign = async (
    client: Socket,
    input: string,
    optionInfo?: SnpPrimerDesignInfo,
  ) => {
    if (this.chatsService.primerDesigns.has(this.conversationUUID)) {
      const previewFunctionCall = {
        name: AgentFunctionByType[AgentType.PRIMER_DESIGN],
        arguments: `{"description":"${input}"}`,
      };
      const agent: BaseAgents = this.functionToAgents[previewFunctionCall.name];
      const agentMessage = await this.chatsService.createMessage(
        '',
        '',
        this.conversationUUID,
        Role.Assistant,
        agent.agent.agentType,
      );
      this.chatsService.emitToCliet(client, TEXT_ANSWER_CREATE, {
        success: true,
        data: agentMessage,
        finishReason: null,
        conversationUUID: this.conversationUUID,
      });
      for await (const msg of this._callFunction(
        previewFunctionCall,
        input,
        optionInfo,
      )) {
        if (msg.role === Role.Assistant) {
          if (msg?.optionInfo) {
            agentMessage.optionInfo = msg.optionInfo;
          }
          agentMessage.content += msg.content;
          this.chatsService.emitToCliet(client, TEXT_ANSWER_GENERATING, {
            success: true,
            data: agentMessage,
            finishReason: null,
            conversationUUID: this.conversationUUID,
          });
          await this.dataSource.manager.update(
            MessageEntity,
            agentMessage.id,
            agentMessage,
          );
          this.chatsService.emitToCliet(client, TEXT_ANSWER_GENERATING, {
            success: true,
            data: agentMessage,
            finishReason: 'stop',
            conversationUUID: this.conversationUUID,
          });
          if (msg?.optionInfo.state == 'stop') {
            await this.handlePrimerDesignStop(client, msg);
          }
        }
      }
    }
  };

  async handleProtocolDesignStop(
    client: Socket,
    responses: {
      optionInfo: ProtocolDesignInfo;
      content: string;
      role: Role;
    },
  ) {
    if (responses.optionInfo?.stage > 1) {
      const responsesStr = JSON.stringify(responses);
      if (responses.optionInfo?.data?.json_file.length > 0) {
        const history = await this.chatsService.findAllMessages(
          this.conversationUUID,
        );
        const currentStep = this.planner.getCurrentStep();
        this.updateStepResult(currentStep.stepName, responsesStr);
        await this.showSummary(client, history);
        this.planner.goNextStep();
      } else {
        const achieved = await this.checkAchieved(client, responsesStr);
        if (!achieved) {
          this.logger.debug('=====ProtocolDesignStop=checkAchieved=fail==');
        }
      }
      await this.chatsService.createMessage(
        INITIATIVE_START_CODE_EXECUTION,
        'We have now completed the Protocol Design step. Do you want to start the code execution process now?',
        this.conversationUUID,
        Role.Assistant,
      );
      client.emit(UPDATE_MESSAGE, {
        conversationUUID: this.conversationUUID,
      });
    } else {
      // search step finished, start normal chat
      this.logger.debug(
        'sendToProtocolDesign==>Protocol Design end, start normal chat',
      );
    }
  }

  sendToProtocolDesign = async (
    client: Socket,
    input: string,
    optionInfo?: any,
  ) => {
    if (this.chatsService.protocolDesigns.has(this.conversationUUID)) {
      const previewFunctionCall = {
        name: AgentFunctionByType[AgentType.PROTOCOL_DESIGN],
        arguments: `{"description":"${input}"}`,
      };
      const agent: BaseAgents = this.functionToAgents[previewFunctionCall.name];
      const agentMessage = await this.chatsService.createMessage(
        '',
        '',
        this.conversationUUID,
        Role.Assistant,
        agent.agent.agentType,
      );
      this.chatsService.emitToCliet(client, TEXT_ANSWER_CREATE, {
        success: true,
        data: agentMessage,
        finishReason: null,
        conversationUUID: this.conversationUUID,
      });
      for await (const msg of this._callFunction(
        previewFunctionCall,
        input,
        optionInfo,
      )) {
        if (msg.role == Role.Assistant) {
          if (msg?.optionInfo) {
            agentMessage.optionInfo = msg.optionInfo;
          }
          agentMessage.content += msg.content;
          this.chatsService.emitToCliet(client, TEXT_ANSWER_GENERATING, {
            success: true,
            data: agentMessage,
            finishReason: null,
            conversationUUID: this.conversationUUID,
          });
          await this.dataSource.manager.update(
            MessageEntity,
            agentMessage.id,
            agentMessage,
          );
          this.chatsService.emitToCliet(client, TEXT_ANSWER_GENERATING, {
            success: true,
            data: agentMessage,
            finishReason: 'stop',
            conversationUUID: this.conversationUUID,
          });
          if (msg?.optionInfo.state == 'stop') {
            await this.handleProtocolDesignStop(client, msg);
          }
        }
      }
    }
  };

  private async handleStreamMessage(
    client: Socket,
    assistantMessage: MessageEntity,
    events:
      | AsyncIterable<ChatCompletions>
      | AsyncIterable<OpenAI.ChatCompletionChunk>,
  ) {
    for await (const event of events) {
      for (const choice of event.choices) {
        let finishReason = '';
        if ('finishReason' in choice) {
          finishReason = choice.finishReason;
        } else {
          finishReason = choice.finish_reason;
        }
        let delta = null;
        if (`delta` in choice) {
          delta = choice.delta;
        }
        if (delta.content) {
          delta = delta.content;
          assistantMessage.content = assistantMessage.content + delta;
          // update database
          await this.dataSource.manager.update(
            MessageEntity,
            assistantMessage.id,
            assistantMessage,
          );
          this.chatsService.emitToCliet(client, TEXT_ANSWER_GENERATING, {
            success: true,
            data: assistantMessage,
            finishReason: finishReason,
            conversationUUID: this.conversationUUID,
          });
        }
        if (finishReason == 'stop') {
          this.chatsService.emitToCliet(client, TEXT_ANSWER_GENERATING, {
            success: true,
            data: assistantMessage,
            finishReason: finishReason,
            conversationUUID: this.conversationUUID,
          });
          // update database
          await this.dataSource.manager.update(
            MessageEntity,
            assistantMessage.id,
            assistantMessage,
          );
        }
      }
    }
  }

  private async checkIfRestart(input: string) {
    let content = `
    # Purpose: Translate user input into English first, then judge "UserInput" means 'restart', 'reboot', 'relaunch', 'reload', 'renew', 'reset', 'regenerate', 'stop', 'finish', 'abort', 'end', 'cancel', 'terminate'. if so, output 'restart' should be true. Especially if UserInput 'start' something, output 'restart' should be false;
    # UserInput:${input};
    # Output: Please  output a string of a JSON object({"restart":true,"thought":""},Don't start with \`\`\`json and end width \`\`\`),including the following key:
    - restart // Boolean type, Whether to achieve the purpose
    - thought // String type, thinking process
    #Rule:
    * You only need to judge directly from what the user says, Don't reason yourself, just understand it directly from the literal meaning
    * Input is what the user says, out is your output
    `;
    // # Examples:
    // - input:start the conversation again,output:restart should be true
    // - input:start the dialogue again,output:restart should be true
    // - input:back to the beginning,output:restart should be true
    // - input:start Primer Design,output:restart should be false
    // - input:all params use default value,output:{"restart":false,"thought":""};

    // if (input === 'restart') {
    //   return true;
    // } else {
    //   return false;
    // }
    content = removeLeadingSpaces(content);
    const startTime = performance.now();
    let events = null;
    if (this.llm_type == LLMType.Azure) {
      const messages = [{ role: Role.User, content }];
      events = await this.openAIClient.getChatCompletions(
        process.env.OPENAI_API_ENGINE,
        messages,
      );
      this.logger.log(`this.openAIClient.events: ${JSON.stringify(events)}`);
      this.updateTokens(events.usage);
    } else {
      const params: OpenAI.ChatCompletionCreateParams = {
        model: process.env.OPENAI_API_ENGINE,
        messages: [{ role: Role.User, content }],
        stream: false,
      };
      events = await this.GPTClient.chat.completions.create(params);
    }
    const endTime = performance.now();
    const executionTime = ((endTime - startTime) / 1000).toFixed(2);
    this.logger.warn(`Execution checkIfRestart time: ${executionTime} seconds`);
    let json = { restart: false };
    if (events.choices && events.choices.length > 0) {
      const choice = events.choices[0];
      if (choice.message) {
        this.logger.debug(`checkIfRestart content: ${choice.message.content}`);
        try {
          json = JSON.parse(choice.message.content);
        } catch {
          json.restart = choice.message.content.includes('true');
        }
      }
    }
    return json.restart;
  }

  private async _handleChatCompletions(
    client: Socket,
    text: string,
    assistantMessage: MessageEntity,
    messages: ChatMessage[],
    events:
      | AsyncIterable<ChatCompletions>
      | AsyncIterable<OpenAI.ChatCompletionChunk>,
  ) {
    const functionCall: FunctionCall = {
      name: '',
      arguments: '',
    };
    for await (const event of events) {
      for (const choice of event.choices) {
        let finishReason = '';
        if ('finishReason' in choice) {
          finishReason = choice.finishReason;
        } else {
          finishReason = choice.finish_reason;
        }
        // text chat
        let fc = null;
        if ('functionCall' in choice.delta) {
          fc = choice.delta.functionCall;
        } else {
          fc =
            (choice.delta as any).function_call ??
            (choice.delta as any).tool_calls;
        }
        if (choice.delta?.content) {
          const delta = choice.delta.content;
          // this.logger.log(`${assistantMessage.id} Chatbot: ${delta}`);
          assistantMessage.content = assistantMessage.content + delta;
          // update database
          await this.dataSource.manager.update(
            MessageEntity,
            assistantMessage.id,
            assistantMessage,
          );
          this.chatsService.emitToCliet(client, TEXT_ANSWER_GENERATING, {
            success: true,
            data: assistantMessage,
            finishReason: finishReason,
            conversationUUID: this.conversationUUID,
          });
        } else if (fc) {
          // function call
          if (fc.name) {
            functionCall.name += fc.name;
          }
          if (fc.arguments) {
            functionCall.arguments += fc.arguments;
          }
        }

        // “stop”, “length”, “content_filter”, “function_call” are finished state
        switch (finishReason) {
          case 'stop': {
            this.logger.log(
              `choice.finishReason stop :${assistantMessage.id} Chatbot: ${assistantMessage.content}`,
            );
            this.chatsService.emitToCliet(client, TEXT_ANSWER_GENERATING, {
              success: true,
              data: assistantMessage,
              finishReason: finishReason,
              conversationUUID: this.conversationUUID,
            });
            break;
          }
          case 'function_call': {
            this.logger.log(
              `choice.finishReason function_call :${assistantMessage.id} Chatbot: ${functionCall.name}`,
            );
            const agent: BaseAgents = this.functionToAgents[functionCall.name];
            // assistantMessage.agentType = agent.agent.agentType;
            // prepare
            const { content: response } = this._prepareCallFunction(
              functionCall,
              text,
            );

            this.logger.debug('function_call', functionCall);
            this.logger.log(
              `${functionCall.name} FunctionCall response: ${response}`,
            );
            // adding assistant response to messages
            let assistantResponse = null;
            if (this.llm_type === LLMType.Azure) {
              assistantResponse = {
                role: Role.Assistant,
                functionCall: { ...functionCall },
                content: undefined,
              };
            } else {
              assistantResponse = {
                role: Role.Assistant,
                function_call: { ...functionCall },
                content: undefined,
              };
            }
            messages.push(assistantResponse);
            // adding function response to messages
            const functionResponse: ChatMessage = {
              role: Role.Function,
              name: functionCall.name,
              content: response,
            };
            messages.push(functionResponse);

            // clear function call
            const previewFunctionCall = { ...functionCall };

            // choose LLM according to llm_type
            let events = null;
            if (this.llm_type === LLMType.Azure) {
              events = await this.openAIClient.listChatCompletions(
                process.env.OPENAI_API_ENGINE,
                messages,
              );
              this.logger.log(
                `this.openAIClient.events: ${JSON.stringify(events)}`,
              );
              this.updateTokens(events.usage);
            } else if (this.llm_type === LLMType.OpenAI) {
              const chatCompletionMessages: OpenAI.ChatCompletionMessageParam[] =
                messages.map((message) => {
                  if ('name' in message) {
                    return message;
                  } else {
                    return { ...message, name: '' }; // or some other default value
                  }
                }) as OpenAI.ChatCompletionMessageParam[];
              const params: OpenAI.ChatCompletionCreateParams = {
                model: process.env.OPENAI_API_ENGINE,
                messages: chatCompletionMessages,
                stream: true,
              };
              this.logger.log('OpenAI ChatCompletionCreateParams', params);
              events = await this.GPTClient.chat.completions.create(params);
            } else {
              events = null;
            }

            await this._handleChatCompletions(
              client,
              text,
              assistantMessage,
              messages,
              events,
            );

            this.logger.log('===============');
            this.logger.log(functionCall.name);
            this.logger.log('===============');
            // create agent message
            const agentMessage = await this.chatsService.createMessage(
              '',
              '',
              this.conversationUUID,
              Role.Assistant,
              agent.agent.agentType,
            );
            this.chatsService.emitToCliet(client, TEXT_ANSWER_CREATE, {
              success: true,
              data: agentMessage,
              finishReason: null,
              conversationUUID: this.conversationUUID,
            });

            for await (const msg of this._callFunction(
              previewFunctionCall,
              text,
              // {},
            )) {
              this.logger.debug('call function msg', msg);
              if (msg.role === Role.Assistant) {
                // only assistant message return to AI context
                if (msg?.optionInfo) {
                  agentMessage.optionInfo = msg.optionInfo;
                }
                agentMessage.content += msg.content;
                this.chatsService.emitToCliet(client, TEXT_ANSWER_GENERATING, {
                  success: true,
                  data: agentMessage,
                  finishReason: null,
                  conversationUUID: this.conversationUUID,
                });
                //execute_ot2_protocol
                if (
                  previewFunctionCall.name ===
                  AgentFunctionByType[AgentType.CODE_EXECUTION]
                ) {
                  const flag = this.cacheService.getObject<number>(
                    `code-execute:${this.conversationUUID}`,
                  );
                  if ((await flag) === 1) {
                    const uploadPcrMessage =
                      await this.chatsService.createMessage(
                        UPLOAD_PRIMER_EXCEL,
                        'please upload primer excel!',
                        this.conversationUUID,
                        Role.Assistant,
                      );

                    client.emit(UPLOAD_PRIMER_EXCEL, {
                      success: true,
                      data: uploadPcrMessage,
                      finishReason: 'stop',
                    });
                  }
                }
              } else if (msg.role === Role.None) {
                if (msg.type === MiddleMessageType.OTAnalysis) {
                  agentMessage.protocolAnalysis = msg.content as any;
                  // update database
                  await this.dataSource.manager.update(
                    MessageEntity,
                    agentMessage.id,
                    agentMessage,
                  );
                  this.chatsService.emitToCliet(
                    client,
                    TEXT_ANSWER_GENERATING,
                    {
                      success: true,
                      data: agentMessage,
                      finishReason: null,
                      conversationUUID: this.conversationUUID,
                    },
                  );
                } else if (msg.type === MiddleMessageType.OTCurrentCommand) {
                  agentMessage.currentCommandIndex = msg.content as any;
                  this.chatsService.emitToCliet(
                    client,
                    TEXT_ANSWER_GENERATING,
                    {
                      success: true,
                      data: agentMessage,
                      finishReason: null,
                      conversationUUID: this.conversationUUID,
                    },
                  );
                } else if (msg.type === MiddleMessageType.JetsonFaultCheck) {
                  agentMessage.faultCheckResult = msg.content as any;
                  this.chatsService.emitToCliet(
                    client,
                    TEXT_ANSWER_GENERATING,
                    {
                      success: true,
                      data: agentMessage,
                      finishReason: null,
                      conversationUUID: this.conversationUUID,
                    },
                  );
                }
              } else {
                this.logger.warn(`unknown msg: ${JSON.stringify(msg)}`);
              }
            }
            await this.dataSource.manager.update(
              MessageEntity,
              agentMessage.id,
              agentMessage,
            );
            this.chatsService.emitToCliet(client, TEXT_ANSWER_GENERATING, {
              success: true,
              data: agentMessage,
              finishReason: 'stop',
              conversationUUID: this.conversationUUID,
            });
            messages.push(agentMessage);

            break;
          }
        }
      }
    }
  }

  private _prepareCallFunction(functionCall: FunctionCall, text: string) {
    if (this.functionToAgents[functionCall.name]) {
      const agent: BaseAgents = this.functionToAgents[functionCall.name];
      return {
        role: Role.Assistant,
        content: `Prompt user, You only need to provide users with a brief explanation of why they should use this agent，\n
          then you will arrange the ${agent.agent.name} \n
          (prompt: ${agent.agent.name}'s description is "${agent.agent.desc}"). \n
          and then tell him the agent is working for him, please wait. \n
          The result will send to user later \n
          Note! The term "Agent" here does not refer to a tool but to an agent or representative.
          Please keep your responses as concise as possible`,
      };
    } else {
      this.logger.log(
        `No matched Agent: ${functionCall.name} response success.`,
      );
      return { role: Role.Assistant, content: `execute successfully` };
    }
  }

  private async *_callFunction(
    functionCall: FunctionCall,
    text: string,
    optionInfo?: any,
  ) {
    if (this.functionToAgents[functionCall.name]) {
      const agent: BaseAgents = this.functionToAgents[functionCall.name];

      this.logger.log(
        `Call Agent: "${agent.agent.name}", args: ${
          functionCall.arguments
        }, optionInfo: ${JSON.stringify(optionInfo)}`,
      );
      this.logger.debug('functionCall.arguments', functionCall.arguments);
      let badJson = false;
      let description = '';
      try {
        const ret = JSON.parse(functionCall.arguments.replaceAll(`\n`, ''));
        description = ret.description;
      } catch (e) {
        badJson = true;
        this.logger.error(`${e.message}`, e.stack);
      }
      if (!text || badJson) {
        this.logger.debug('_callFunction ==> description', description);
        this.logger.debug('_callFunction ==> text', text);
        this.logger.debug('_callFunction ==> badJson', badJson);
        yield {
          role: Role.Assistant,
          content: `Argument parse error，please try again.`,
        };
        return;
      }

      const stepInfo = Object.keys(this.cvs.stepsResult)
        .map((stepName) => {
          return `${stepName}: ${JSON.stringify(
            this.cvs.stepsResult[stepName],
          )}`;
        })
        .join('\n');
      const currentStep = this.planner.getCurrentStep();
      const isDeterminingDNA =
        currentStep.stepName == this.planner.steps[1].stepName;
      // add text to description
      if (!isDeterminingDNA) {
        text = `
          ${stepInfo}
          ${text}
        `;
      }

      yield* agent.send({
        userInput: text,
        description,
        optionInfo,
      });
    } else {
      this.logger.log(
        `No matched Agent: ${functionCall.name} response success.`,
      );
      yield { role: Role.Assistant, content: `execute successfully, no error` };
    }
  }
}
