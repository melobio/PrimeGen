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
import { DataSource } from 'typeorm';
import { ExperimentEntity } from '../../experiments/entities/experiment.entity';
import {
  AgentFunctionByType,
  AgentFunctions,
  AgentType,
  LLMFields,
  LLMType,
  NodeTypes,
} from '@xpcr/common';
import { NCBIFunctionType, NCBIResultContent } from '../agents/ext/ncbi-search';
import { LlmsEntity } from '../../tools/llms/entities/llms.entity';
import { Agents } from '../../tools/agents/entities/agents.entity';
import { FaultAgents } from './fault-agents';
import { CodeExecutionAgents } from './code-execution-agents';
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
} from '../conversation.constant';
import { InternetSearchAgents } from './internet-search-agents';
import { SequenceSearchAgents } from './sequence-search-agents';
import { ConversationService } from '../conversation.service';
import { PlanStep, Planner } from './ext/planner';
import { ConversationEntity } from '../entities/conversation.entity';
import { PrimerDesignAgents } from './primer-design-agents';
import { ProtocolDesignAgents } from './protocol-design-agents';
import { CacheService } from 'src/cache/cache.service';
interface jasonType {
  output: string;
  achieved: boolean;
  thought: string;
  objective?: string;
  wantTips?: string;
}
function removeLeadingSpaces(str) {
  // 使用正则表达式替换每一行开头的空格
  return str.replace(/^\s+/gm, '');
}

function removeLineBreaks(str) {
  // 使用正则表达式替换换行符
  return str.replace(/\n/g, '');
}
// 正则判断search 是否达到目的(输出序列文件或提供基因信息)
function checkStringForFileOrGeneInfo(inputString) {
  console.log(`checkStringForFileOrGeneInfo inputString: ${inputString}`);
  // 正则表达式用于匹配文件名
  const fileRegex = /\.(fasta|fastq|fna|csv)\b/;

  // 正则表达式用于匹配基因信息字段
  const geneInfoRegex1 =
    /\b(?=.*chrom)(?=.*start)(?=.*end)(?=.*gene)(?=.*type)(?=.*HGNC)\b/;
  const geneInfoRegex2 =
    /\b(?=.*Gene name)(?=.*Chromosome)(?=.*Gene start)(?=.*Gene end)(?=.*Start position in gene)(?=.*End position in gene)\b/;
  // 检查输入字符串是否匹配文件名或基因信息字段
  if (
    fileRegex.test(inputString) ||
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
  nextMaybeGoNextStep = false;
  planner: Planner;
  openAIClient: OpenAIClient;
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

  // 更新步骤为StepName的结果
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

    // 读取check-purpose模板内容
    this.checkPurpose = fs
      .readFileSync(
        path.join(process.cwd(), 'assets', 'plans', 'check-purpose.txt'),
      )
      .toString('utf-8');
    this.logger.debug(`checkPurpose: ${String(this.checkPurpose)}`);

    // 读取next模板内容
    this.next = fs
      .readFileSync(path.join(process.cwd(), 'assets', 'plans', 'next.txt'))
      .toString('utf-8');
    this.logger.debug(`next: ${String(this.next)}`);

    // 读取summary模板内容
    this.summary = fs
      .readFileSync(path.join(process.cwd(), 'assets', 'plans', 'summary.txt'))
      .toString('utf-8');
    this.logger.debug(`summary: ${String(this.summary)}`);

    // 读取PCR实验的预设数据,并初始化planner,并设置cvs的currentStep
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
    this.openAIClient = new OpenAIClient(
      apiBase,
      new AzureKeyCredential(apiKey),
    );
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
            newAgent = new CodeExecutionAgents(
              agent,
              this.dataSource,
              this.cacheService,
              this.chatsService,
            );
            break;
          case AgentType.INTERNET_SEARCH:
            newAgent = new InternetSearchAgents(agent, this.dataSource);
            break;
          case AgentType.SEQUENCE_SEARCH:
            newAgent = new SequenceSearchAgents(
              agent,
              this.dataSource,
              this.chatsService,
            );
            break;
          case AgentType.PRIMER_DESIGN:
            newAgent = new PrimerDesignAgents(agent, this.dataSource);
            break;
          case AgentType.PROTOCOL_DESIGN:
            newAgent = new ProtocolDesignAgents(agent, this.dataSource);
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
          // 一个字符串传递Agent需要做的事情，TODO 应对每种Agent增加特定参数，更准确
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

  // 完成了STEP的第一步，需要提示用户下一步的步骤
  async showSummary(client: Socket, history: MessageEntity[]) {
    const nextStep = this.planner.getNextStep();

    const assistantMessage = await this.chatsService.createMessage(
      '',
      '',
      this.conversationUUID,
      Role.Assistant,
    );
    // 发送给客户端
    client.emit(TEXT_ANSWER_CREATE, {
      success: true,
      data: assistantMessage,
      finishReason: null,
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
    const messages = [
      { role: Role.System, content: system },
      { role: Role.User, content },
    ];
    const startTime = performance.now();
    this.logger.debug(`summary messages:${JSON.stringify(messages)}`);
    const events = await this.openAIClient.listChatCompletions(
      process.env.OPENAI_API_ENGINE,
      messages,
    );
    await this.handleStreamMessage(client, assistantMessage, events);

    const endTime = performance.now();
    const executionTime = ((endTime - startTime) / 1000).toFixed(2);
    this.logger.warn(`Execution summary time: ${executionTime} seconds`);

    // if (events.choices && events.choices.length > 0) {
    //   const choice = events.choices[0];

    //   if (choice.message) {
    //     assistantMessage.content = choice.message.content;
    //     this.logger.warn(`Summary message: ${choice.message.content}`);
    //   }
    // }

    ////////////////////// auto next
    // let canGoNext = false;
    // let goNextText = '';
    // content = `
    //   #目标: 帮助用户判断是否能进行下一步实验
    //   #思考: 分析实验的SOP和提供的信息，判断是否具备了进行下一步实验的条件，如果具备了条件，那么将进行下一步。
    //   #输出：请输出一个JSON对象的字符串，包括以下字段:
    //   - goNextText // String类型，用一个句子描述，描述的内容类似 "进行xxx"，其中xxx表示下一步具体的内容
    //   - canGoNext // Boolean类型，是否可以进行下一步
    //   - thought // String类型，思考过程
    // `;
    // content = removeLeadingSpaces(content);
    //
    // this.logger.debug(`auto next content: ${content}`);
    //
    // chatHistory.push({
    //   role: Role.Assistant,
    //   content: assistantMessage.content,
    // });
    //
    // messages = [
    //   { role: Role.System, content: system },
    //   ...chatHistory,
    //   { role: Role.User, content },
    // ];
    //
    // events = await this.openAIClient.getChatCompletions(
    //   process.env.OPENAI_API_ENGINE,
    //   messages,
    // );
    //
    // if (events.choices && events.choices.length > 0) {
    //   const choice = events.choices[0];
    //
    //   if (choice.message) {
    //     console.log(messages);
    //     // assistantMessage.content = choice.message.content;
    //     const text = choice.message.content.replace(/```json|```/g, '');
    //     this.logger.warn(`Auto next message: ${text}`);
    //     try {
    //       const json = JSON.parse(text);
    //       canGoNext = json.canGoNext;
    //       goNextText = json.goNextText;
    //     } catch (e) {
    //       this.logger.warn(`Auto next JSON parse error: ${e.message}`);
    //     }
    //   }
    // }

    //////////////////////

    /////////////////// auto next
    // if (canGoNext) {
    //   history = [...history, assistantMessage];
    //   assistantMessage = await this.chatsService.createMessage(
    //     '',
    //     '',
    //     this.conversationUUID,
    //     Role.Assistant,
    //   );
    //   // 发送给客户端
    //   client.emit(TEXT_ANSWER_CREATE, {
    //     success: true,
    //     data: assistantMessage,
    //     finishReason: null,
    //   });
    //
    //   await this.sendText(
    //     client,
    //     { content: goNextText, role: Role.User } as MessageEntity,
    //     assistantMessage,
    //     history,
    //   );
    // }
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
    const chatHistory = history.map<ChatMessage>((item) => {
      return { role: item.role, content: item.content };
    });
    let json: jasonType = {
      output: '',
      achieved: false,
      thought: '',
      objective: '',
      wantTips: '',
    };
    const cvs = history
      .map((item) => `${item.role}: ${item.content}`)
      .slice(1)
      .join('\n');
    const userFirstMsg = chatHistory.find((item) => item.role == Role.User);
    // 就提取具体的实验目的(具体的实验类型和实验对象)是什么/具体的文件列表是什么/具体的Protocol脚本等等
    let content = `
    #Purpose: Assistant needs to quickly get "${
      currentStep.purpose
    }" from the dialogue
    #Think: Analyze the dialogue, whether the assistant can get the purpose from the dialogue
    #Output：Please quickly output a string of a JSON object(like {"output":"","achieved":"","thought":"","objective":"","wantTips":""}),including the following key:
    - output // String type, the purpose summarized from the dialogue, the purpose of the dialogue needs to include the actual results.
    - achieved // Boolean type,Whether to achieve the purpose
    - thought // String type, thinking process
    - wantTips // String type, when the purpose is not achieved, the assistant needs to provide a prompt to achieve the purpose of the user to achieve (please explain ...  So I can better help you complete ${
      currentStep.purpose
    });
    ${
      isConfirmExperObj
        ? `- objective // String type, experimental purpose(words composed of experimental objectives and experiment types);(The first sentence that the user says is${userFirstMsg.content}Please determine the language output of the following fields according to the first sentence input of the user, default English output);`
        : ''
    }
    # Rule
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

    const messages = [
      { role: Role.System, content: system },
      { role: Role.User, content },
    ];
    const startTime = performance.now();
    this.logger.warn(`checkAchieved messages:${JSON.stringify(messages)}`);
    const events = await this.openAIClient.getChatCompletions(
      process.env.OPENAI_API_ENGINE,
      messages,
    );
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
      await this.showSummary(client, history);
      this.planner.goNextStep();
      if (isConfirmExperObj && json.objective) {
        this.logger.debug('sendTextToSearch====>0');
        await this.sendTextToSearch(client, input, chatHistory);
        // 更新会话名称
        this.cvs.name = json.objective;
        this.dataSource.manager
          .update(ConversationEntity, this.cvs.id, this.cvs)
          .then();
        // 更新实验名称
        const experiment = await this.dataSource.manager.findOne(
          ExperimentEntity,
          {
            where: {
              conversationUUID: this.cvs.uuid,
            },
          },
        );
        experiment.name = `${json.objective} With OT-2`;
        this.dataSource.manager
          .update(ExperimentEntity, experiment.id, experiment)
          .then();
      } else {
        this.logger.debug('TEXT_ANSWER_DONE 1');
        client.emit(TEXT_ANSWER_DONE, {
          conversationUUID: this.conversationUUID,
        });
      }
    } else {
      // if (json.wantTips) {
      //   this.handleTipsAfterCheckAchieved(client, json.wantTips);
      // }
      this.logger.debug('TEXT_ANSWER_DONE 2');
      client.emit(TEXT_ANSWER_DONE, {
        conversationUUID: this.conversationUUID,
      });
    }
  }

  async handleNcbiSearchStop(
    client: Socket,
    responses: NCBIResultContent['responses'],
  ) {
    const history = await this.chatsService.findAllMessages(
      this.conversationUUID,
    );
    const currentStep = this.planner.getCurrentStep();
    if (responses.stage > 2) {
      // search阶段结束，开始总结
      this.logger.debug('sendTextToSearch ==>search阶段结束，开始反思');
      this.updateStepResult(currentStep.stepName, JSON.stringify(responses));
      await this.showSummary(client, history);
      this.planner.goNextStep();
      // 下一步引物设计，创建开始引物设计的提示气泡
      this.createPrimerDesignTips(client);
    } else {
      // search阶段结束，开始闲聊
      this.logger.debug('sendTextToSearch ==>search阶段结束，开始闲聊');
    }
  }
  async createPrimerDesignTips(client: Socket) {
    // 创建start PrimerDesign提示消息
    const tipsMessage = await this.chatsService.createMessage(
      'INITIATIVE_START_PRIMER_DESIGN',
      'We have now completed the Determing DNA Sequence step. Do you want to start the primer design process now?',
      this.conversationUUID,
      Role.Assistant,
    );

    client.emit(TEXT_ANSWER_CREATE, {
      success: true,
      data: tipsMessage,
      finishReason: null,
    });
  }

  async handleTipsAfterCheckAchieved(client: Socket, tipsText: string) {
    // 创建一个Assistant消息, 提示缺乏什么达成目标
    const assistantMessage = await this.chatsService.createMessage(
      '',
      tipsText,
      this.conversationUUID,
      Role.Assistant,
    );
    client.emit(TEXT_ANSWER_CREATE, {
      success: true,
      data: assistantMessage,
      finishReason: null,
    });
    client.emit(TEXT_ANSWER_GENERATING, {
      success: true,
      data: assistantMessage,
      finishReason: 'stop',
      conversationUUID: this.conversationUUID,
    });
  }

  async sendText(
    client: Socket,
    userMessageContent: string,
    history: MessageEntity[],
  ) {
    const content = removeLeadingSpaces(userMessageContent);
    // if (this.nextMaybeGoNextStep) {
    //   this.nextMaybeGoNextStep = false;
    //   // await this.getCurrentStepActualResult(history);
    //   if (await this.goNextStep(client, content, assistantMessage, history)) {
    //     return;
    //   }
    // }

    // content = `
    //   # Q&A Examples
    //   ${this.examples[this.planner.getCurrentStep().stepName].replace(
    //     '{{NEXT_STEP}}',
    //     this.planner.getNextStep().stepName,
    //   )}
    //
    //   # User Input: ${content}, (Reply with reference to the Q&A Examples, Do not mention the word "example" in the answer.)
    // `;

    this.logger.log(`content: ${content}`);

    const chatHistory = history.map<ChatMessage>((item) => {
      return { role: item.role, content: item.content };
    });

    const currentStepIndex = this.planner.steps.indexOf(
      this.planner.getCurrentStep(),
    );
    const doneSteps =
      currentStepIndex !== 0
        ? this.planner.steps.slice(0, currentStepIndex + 1)
        : [];
    this.initMessages.forEach((item, index) => {
      if (item.role == Role.System) {
        this.initMessages[index].content = this.initMessages[
          index
        ].content.replace('{{DONE_STEPS}}', JSON.stringify(doneSteps));
      }
    });
    const messages = [
      ...this.initMessages,
      ..._.slice(chatHistory, -50), // TODO 应该计算token能覆盖到前面几条记录，而不是直接取最后50条
      {
        role: Role.User,
        content,
      },
    ];
    this.logger.log(`messages: ${JSON.stringify(messages)}`);
    const functions = this._getFunctions();
    this.logger.log(`functions: ${JSON.stringify(functions)}`);
    const currentStep = this.planner.getCurrentStep();
    const isDeterminingDNA =
      currentStep.stepName == this.planner.steps[1].stepName;
    const isPrimerDesign =
      currentStep.stepName == this.planner.steps[2].stepName;
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
    // 确定`基因步骤`有什么agent就直接调用agent，绕过planner
    if (isDeterminingDNA && searchAgentMessages.length > 0) {
      await this.sendTextToSearch(client, content, messages);
    } else if (isPrimerDesign && primerAgentMessages.length > 0) {
      await this.sendTextToPrimerDesign(client, content, messages);
    } else {
      this.logger.debug('=normal_talk=======>');
      const assistantMessage = await this.chatsService.createMessage(
        '',
        '',
        this.conversationUUID,
        Role.Assistant,
      );
      client.emit(TEXT_ANSWER_CREATE, {
        success: true,
        data: assistantMessage,
        finishReason: null,
      });
      const events = await this.openAIClient.listChatCompletions(
        process.env.OPENAI_API_ENGINE,
        messages,
        functions.length > 0
          ? {
              functionCall: 'auto',
              functions: functions,
            }
          : null,
      );

      // try catch 是处理openai的关键词
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
        // 更新数据库
        await this.dataSource.manager.update(
          MessageEntity,
          assistantMessage.id,
          assistantMessage,
        ); // 发送给客户端
        client.emit(TEXT_ANSWER_GENERATING, {
          success: true,
          data: assistantMessage,
          finishReason: null,
          conversationUUID: this.conversationUUID,
        });

        // 发送给客户端
        client.emit(TEXT_ANSWER_GENERATING, {
          success: true,
          data: assistantMessage,
          finishReason: 'stop',
          conversationUUID: this.conversationUUID,
        });
        throw e;
      }

      // 判断是否应该进入下一步，如果需要，提示用户
      // this.logger.debug(`ask whether go next step`);
      // await this.askWhetherGoNextStep(client, [
      //   ...history,
      //   userMessage,
      //   assistantMessage,
      // ]);

      // 判断用户是否已同意进入到下一步
      // const shouldGoNext = await this.checkGoNext(
      //   client,
      //   userMessage.content,
      //   history,
      // );
      // this.logger.debug(`should go next step: ${shouldGoNext}`);
      // if (shouldGoNext) {
      //   this.planner.goNextStep();
      //   await this.startCurrentStep(client, history);
      // }

      // 基于历史对话/agent消息，判断是否达到目标
      await this.checkAchieved(client, content);
    }
  }

  // 直接调用agent，绕过planner
  sendTextToSearch = async (
    client: Socket,
    input: string,
    messages: ChatMessage[],
  ) => {
    const previewFunctionCall = {
      name: AgentFunctionByType[AgentType.SEQUENCE_SEARCH],
      arguments: `{"description":"${input}"}`,
    };
    const agent: BaseAgents = this.functionToAgents[previewFunctionCall.name];
    const history = await this.chatsService.findAllMessages(
      this.conversationUUID,
    );
    let stage = 1;
    let search_type: NCBIFunctionType;
    let protein_mutation_dict = {};
    const searchAgentMessage = history.filter(
      (item) =>
        item.role == Role.Assistant &&
        item.agentType == AgentType.SEQUENCE_SEARCH &&
        item.optionInfo,
    );
    const stageArr = searchAgentMessage
      .map((optItem) => optItem.optionInfo.stage)
      .sort((a, b) => b - a);
    if (stageArr.length > 0) {
      stage = stageArr.pop();
      const optionInfo =
        searchAgentMessage[searchAgentMessage.length - 1].optionInfo;
      search_type = optionInfo?.search_type;
      protein_mutation_dict = optionInfo?.protein_mutation_dict;
    }
    const res = await this.chatsService.handleNcbiSearch({
      option: { stage, search_type, protein_mutation_dict },
      conversationUUID: this.conversationUUID,
    });
    client.emit(TEXT_ANSWER_DONE, {
      conversationUUID: this.conversationUUID,
    });
    if (res.data.state == 'stop' && res.success) {
      await this.handleNcbiSearchStop(client, res.data);
    }
    this.logger.debug('sendTextToSearch====>2');
  };

  sendTextToPrimerDesign = async (
    client: Socket,
    input: string,
    messages: ChatMessage[],
  ) => {
    this.logger.debug('=isPrimerDesign=======>');
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

    client.emit(TEXT_ANSWER_CREATE, {
      success: true,
      data: agentMessage,
      finishReason: null,
    });
    for await (const msg of this._callFunction(previewFunctionCall, input)) {
      if (msg.role === Role.Assistant) {
        // 只有assistant的消息才会返回给AI上下文
        if (msg?.optionInfo) {
          agentMessage.optionInfo = msg.optionInfo;
        }
        agentMessage.content += msg.content;
        client.emit(TEXT_ANSWER_GENERATING, {
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
        client.emit(TEXT_ANSWER_GENERATING, {
          success: true,
          data: agentMessage,
          finishReason: 'stop',
          conversationUUID: this.conversationUUID,
        });

        messages.push(agentMessage);
      }
    }
  };

  private async handleStreamMessage(
    client: Socket,
    assistantMessage: MessageEntity,
    events: AsyncIterable<ChatCompletions>,
  ) {
    for await (const event of events) {
      for (const choice of event.choices) {
        if (choice.delta?.content) {
          const delta = choice.delta.content;
          assistantMessage.content = assistantMessage.content + delta;
          // 更新数据库
          await this.dataSource.manager.update(
            MessageEntity,
            assistantMessage.id,
            assistantMessage,
          );
          // 发送给客户端
          client.emit(TEXT_ANSWER_GENERATING, {
            success: true,
            data: assistantMessage,
            finishReason: choice.finishReason,
            conversationUUID: this.conversationUUID,
          });
        }
        if (choice.finishReason == 'stop') {
          // 发送给客户端
          client.emit(TEXT_ANSWER_GENERATING, {
            success: true,
            data: assistantMessage,
            finishReason: choice.finishReason,
            conversationUUID: this.conversationUUID,
          });

          // 更新数据库
          await this.dataSource.manager.update(
            MessageEntity,
            assistantMessage.id,
            assistantMessage,
          );
        }
      }
    }
  }

  private async _handleChatCompletions(
    client: Socket,
    text: string,
    assistantMessage: MessageEntity,
    messages: ChatMessage[],
    events: AsyncIterable<ChatCompletions>,
  ) {
    const functionCall: FunctionCall = {
      name: '',
      arguments: '',
    };
    for await (const event of events) {
      for (const choice of event.choices) {
        // this.logger.log(choice);

        // text chat
        if (choice.delta?.content) {
          const delta = choice.delta.content;
          // this.logger.log(`${assistantMessage.id} Chatbot: ${delta}`);
          assistantMessage.content = assistantMessage.content + delta;
          // 更新数据库
          await this.dataSource.manager.update(
            MessageEntity,
            assistantMessage.id,
            assistantMessage,
          );
          // 发送给客户端
          client.emit(TEXT_ANSWER_GENERATING, {
            success: true,
            data: assistantMessage,
            finishReason: choice.finishReason,
            conversationUUID: this.conversationUUID,
          });
        } else if (choice.delta?.functionCall) {
          // function call
          const fc = choice.delta.functionCall;
          if (fc.name) {
            functionCall.name += fc.name;
          }
          if (fc.arguments) {
            functionCall.arguments += fc.arguments;
          }
        }

        // “stop”, “length”, “content_filter”, “function_call” 是完成状态
        switch (choice.finishReason) {
          case 'stop': {
            this.logger.log(
              `choice.finishReason stop :${assistantMessage.id} Chatbot: ${assistantMessage.content}`,
            );
            // 发送给客户端
            client.emit(TEXT_ANSWER_GENERATING, {
              success: true,
              data: assistantMessage,
              finishReason: choice.finishReason,
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
            // if (agent.needSummarize) {
            //
            // } else {
            //   // 不需要总结
            //   // 直接发送FunctionCall的结果到Client
            //   await this.dataSource.manager.update(
            //     MessageEntity,
            //     assistantMessage.id,
            //     assistantMessage,
            //   );
            //   client.emit(TEXT_ANSWER_GENERATING, {
            //     success: true,
            //     data: assistantMessage,
            //     finishReason: 'stop',
            //   });
            // }
            // adding assistant response to messages
            const assistantResponse: ChatMessage = {
              role: Role.Assistant,
              functionCall: { ...functionCall },
              content: undefined,
            };
            messages.push(assistantResponse);
            // adding function response to messages
            const functionResponse: ChatMessage = {
              role: Role.Function,
              name: functionCall.name,
              content: response,
            };
            messages.push(functionResponse);
            // messages.push({
            //   role: Role.User,
            //   content:
            //     '基于上面的函数结果进行思考，为什么使用调用这个函数，并对结果进行评估',
            // });

            // clear function call
            const previewFunctionCall = { ...functionCall };
            functionCall.name = '';
            functionCall.arguments = '';

            // const functions = this._getFunctions();
            const events = await this.openAIClient.listChatCompletions(
              process.env.OPENAI_API_ENGINE,
              messages,
              // 避免反复调用function
              // {
              //   functionCall: 'auto',
              //   functions: functions,
              // },
            );
            await this._handleChatCompletions(
              client,
              text,
              assistantMessage,
              messages,
              events,
            );

            // 创建agent消息
            const agentMessage = await this.chatsService.createMessage(
              '',
              '',
              this.conversationUUID,
              Role.Assistant,
              agent.agent.agentType,
            );

            client.emit(TEXT_ANSWER_CREATE, {
              success: true,
              data: agentMessage,
              finishReason: null,
            });

            for await (const msg of this._callFunction(
              previewFunctionCall,
              text,
            )) {
              this.logger.debug('call function msg', msg);
              if (msg.role === Role.Assistant) {
                // 只有assistant的消息才会返回给AI上下文
                if (msg?.optionInfo) {
                  agentMessage.optionInfo = msg.optionInfo;
                }
                agentMessage.content += msg.content;
                client.emit(TEXT_ANSWER_GENERATING, {
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
                  // 更新数据库
                  await this.dataSource.manager.update(
                    MessageEntity,
                    agentMessage.id,
                    agentMessage,
                  );
                  client.emit(TEXT_ANSWER_GENERATING, {
                    success: true,
                    data: agentMessage,
                    finishReason: null,
                    conversationUUID: this.conversationUUID,
                  });
                } else if (msg.type === MiddleMessageType.OTCurrentCommand) {
                  agentMessage.currentCommandIndex = msg.content as any;
                  client.emit(TEXT_ANSWER_GENERATING, {
                    success: true,
                    data: agentMessage,
                    finishReason: null,
                    conversationUUID: this.conversationUUID,
                  });
                } else if (msg.type === MiddleMessageType.JetsonFaultCheck) {
                  agentMessage.faultCheckResult = msg.content as any;
                  // 发送
                  client.emit(TEXT_ANSWER_GENERATING, {
                    success: true,
                    data: agentMessage,
                    finishReason: null,
                    conversationUUID: this.conversationUUID,
                  });
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
            client.emit(TEXT_ANSWER_GENERATING, {
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
      return { role: Role.Assistant, content: `执行完成，没有发现错误` };
    }
  }

  private async *_callFunction(functionCall: FunctionCall, text: string) {
    if (this.functionToAgents[functionCall.name]) {
      const agent: BaseAgents = this.functionToAgents[functionCall.name];

      this.logger.log(
        `Call Agent: "${agent.agent.name}", args: ${functionCall.arguments}`,
      );
      this.logger.debug('functionCall.arguments', functionCall.arguments);
      let badJson = false;
      let description = '';
      try {
        const ret = JSON.parse(functionCall.arguments);
        description = ret.description;
      } catch (e) {
        badJson = true;
        this.logger.error(`${e.message}`, e.stack);
      }
      if (!description || !text || badJson) {
        yield {
          role: Role.Assistant,
          content: `Argument parse error，can not run.`,
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
      // 添加需要的提示到text
      if (!isDeterminingDNA) {
        text = `
          ${stepInfo}
          
          ${text}
        `;
      }

      yield* agent.send(text, description);
    } else {
      this.logger.log(
        `No matched Agent: ${functionCall.name} response success.`,
      );
      yield { role: Role.Assistant, content: `执行完成，没有发现错误` };
    }
  }
}
