import { Injectable } from '@nestjs/common';
import { Agents } from './agents/entities/agents.entity';
import {
  AgentFunctions,
  AgentType,
  DeviceFields,
  DeviceType,
  InputFields,
  InputType,
  LLMFields,
  LLMType,
  NodeFieldType,
  NodeTypes,
  PromptFields,
  PromptType,
  ToolFields,
  ToolType,
} from '@xpcr/common';
import { LlmsEntity } from './llms/entities/llms.entity';
import { DevicesEntity } from './devices/entities/devices.entity';
import { PromptEntity } from './prompts/entities/prompt.entity';
import { InputEntity } from './input/entities/input.entity';
import { ToolEntity } from './tool/entities/tool.entity';

@Injectable()
export class ToolsService {
  async getAgentTemplates() {
    return [
      new Agents(
        AgentType.SEQUENCE_SEARCH,
        `${AgentType.SEQUENCE_SEARCH} Agent`,
        'This sequence search agent function is invoked when a user has a requirement to download bioinformatics data, such as DNA sequences or protein sequences. It is important not to confuse this with general internet searches; this method is based on biological database APIs and web crawling, not on search engines like Google or Wikipedia.',
        [
          {
            name: AgentFunctions.NCBI_SEARCH,
            description:
              'This sequence search agent function is invoked when a user has a requirement to download bioinformatics data, such as DNA sequences or protein sequences. It is important not to confuse this with general internet searches; this method is based on biological database APIs and web crawling, not on search engines like Google or Wikipedia.',
            parameters: {
              type: 'object',
              properties: {
                query: {
                  type: 'string',
                  description: 'A search query which user input.',
                },
              },
            },
            required: ['query'],
          },
        ],
        [],
        [
          {
            type: NodeTypes.LLMs,
            name: 'LLM',
            uuid: '',
            functionName: '',
            functionDesc: '',
            required: true,
          },
        ],
        {
          type: NodeTypes.Agents,
          name: NodeTypes.Agents,
          uuid: '',
          functionName: '',
          functionDesc: '',
          required: false,
        },
        '0.0.1',
      ),
      new Agents(
        AgentType.PRIMER_DESIGN,
        `${AgentType.PRIMER_DESIGN} Agent`,
        '"This is a primer design agent  function. This method is triggered under two conditions: 1) the sequence_search agent has completed its task; 2) the user has confirmed their desire to start primer design for experiments. It communicates with the user to confirm the primer parameters and employs an algorithm to carry out the primer design.',
        [
          {
            name: AgentFunctions.DESIGN_PRIMER,
            description:
              '"This is a primer design agent  function. This method is triggered under two conditions: 1) the sequence_search agent has completed its task; 2) the user has confirmed their desire to start primer design for experiments. It communicates with the user to confirm the primer parameters and employs an algorithm to carry out the primer design.',
            parameters: {
              type: 'object',
              properties: {
                query: {
                  type: 'string',
                  description: 'A search query which user input.',
                },
              },
            },
            required: ['query'],
          },
        ],
        [],
        [
          {
            type: NodeTypes.LLMs,
            name: 'LLM',
            uuid: '',
            functionName: '',
            functionDesc: '',
            required: true,
          },
        ],
        {
          type: NodeTypes.Agents,
          name: NodeTypes.Agents,
          uuid: '',
          functionName: '',
          functionDesc: '',
          required: false,
        },
        '0.0.1',
      ),
      new Agents(
        AgentType.PROTOCOL_DESIGN,
        `${AgentType.PROTOCOL_DESIGN} Agent`,
        'A protocol design toolkit.Useful for when you need to design protocols for experiments.',
        [
          {
            name: AgentFunctions.DESIGN_PROTOCOL,
            description:
              'A protocol design toolkit.Useful for when you need to design protocols for experiments.',
            parameters: {
              type: 'object',
              properties: {
                query: {
                  type: 'string',
                  description: 'A search query which user input.',
                },
              },
            },
            required: ['query'],
          },
        ],
        [],
        [
          {
            type: NodeTypes.LLMs,
            name: 'LLM',
            uuid: '',
            functionName: '',
            functionDesc: '',
            required: true,
          },
        ],
        {
          type: NodeTypes.Agents,
          name: NodeTypes.Agents,
          uuid: '',
          functionName: '',
          functionDesc: '',
          required: false,
        },
        '0.0.1',
      ),
      new Agents(
        AgentType.INTERNET_SEARCH,
        `${AgentType.INTERNET_SEARCH} Agent`,
        `This is a internet search agent function. It is triggered when a user explicitly express for real-time internet information retrieval, such as 'help me retrieve/search  ** information', 'check the weather for me', or 'search the internet for ***'. Note that this search should be distinguished from sequence_search agent; internet_search does not retrieve bioinformatics sequences`,
        [
          {
            name: AgentFunctions.GOOGLE_SEARCH,
            description: `This is a internet search agent function. It is triggered when a user explicitly express for real-time internet information retrieval, such as 'help me retrieve/search  ** information', 'check the weather for me', or 'search the internet for ***'. Note that this search should be distinguished from sequence_search agent; internet_search does not retrieve bioinformatics sequences`,
            parameters: {
              type: 'object',
              properties: {
                query: {
                  type: 'string',
                  description: 'A search query which user input.',
                },
              },
            },
            required: ['query'],
          },
          {
            name: AgentFunctions.WIKIPEDIA_SEARCH,
            description:
              'A tool for interacting with and fetching data from the Wikipedia API.',
            parameters: {
              type: 'object',
              properties: {
                query: {
                  type: 'string',
                  description: 'A search query which user input.',
                },
              },
            },
            required: ['query'],
          },
        ],
        [],
        [
          {
            type: NodeTypes.LLMs,
            name: 'LLM',
            uuid: '',
            functionName: '',
            functionDesc: '',
            required: true,
          },
          {
            type: NodeTypes.Tools,
            name: 'Google Search',
            uuid: '',
            functionName: '',
            functionDesc: '',
            required: true,
          },
          {
            type: NodeTypes.Tools,
            name: 'Wikipedia Search',
            uuid: '',
            functionName: '',
            functionDesc: '',
            required: true,
          },
        ],
        {
          type: NodeTypes.Agents,
          name: NodeTypes.Agents,
          uuid: '',
          functionName: '',
          functionDesc: '',
          required: false,
        },
        '0.0.1',
      ),
      new Agents(
        AgentType.FAULT,
        `${AgentType.FAULT} Agent`,
        'Check whether the OT2 or OT3 machine has encountered a runtime error based on the provided parameters.',
        [
          {
            name: AgentFunctions.CHECK_OT2_STATE,
            description:
              'Check whether the OT2 machine has encountered a runtime error based on the provided parameters.',
            parameters: {
              type: 'object',
              properties: {
                checkType: {
                  type: 'string',
                  description:
                    "The type of the check, 'all', tips' or 'liquid'.",
                },
              },
            },
            required: ['checkType'],
          },
        ],
        [],
        [
          {
            type: NodeTypes.LLMs,
            name: 'LLM',
            uuid: '',
            functionName: '',
            functionDesc: '',
            required: true,
          },
          {
            type: NodeTypes.Devices,
            name: 'OT2',
            uuid: '',
            functionName: '',
            functionDesc: '',
            required: true,
          },
        ],
        {
          type: NodeTypes.Agents,
          name: NodeTypes.Agents,
          uuid: '',
          functionName: '',
          functionDesc: '',
          required: false,
        },
        '0.0.1',
      ),
      new Agents(
        AgentType.CODE_EXECUTION,
        `${AgentType.CODE_EXECUTION} Agent`,
        'OT2 protocol executor. Please input a complete Python script code that can be accurately executed on OT-2',
        [
          {
            name: AgentFunctions.EXECUTE_OT2_PROTOCOL,
            description:
              'Code Execution Agent is a software component designed to execute the protocol on OT2 and return the results.',
            parameters: {
              type: 'object',
              properties: {
                protocol: {
                  type: 'string',
                  description:
                    'The opentrons protocol to be executed, empty when not provided',
                },
              },
            },
            required: ['protocol'],
          },
        ],
        [],
        [
          {
            type: NodeTypes.LLMs,
            name: 'LLM',
            uuid: '',
            functionName: '',
            functionDesc: '',
            required: true,
          },
          {
            type: NodeTypes.Devices,
            name: 'OT2',
            uuid: '',
            functionName: '',
            functionDesc: '',
            required: true,
          },
        ],
        {
          type: NodeTypes.Agents,
          name: NodeTypes.Agents,
          uuid: '',
          functionName: '',
          functionDesc: '',
          required: false,
        },
        '0.0.1',
      ),
    ];
  }

  async getLLMTemplates() {
    // this.logger.log(`open ai key: ${process.env.OPENAI_API_KEY}`);
    // this.logger.log(`open ai base: ${process.env.OPENAI_API_BASE}`);
    // this.logger.log(`open ai engine: ${process.env.OPENAI_API_ENGINE}`);
    return [
      new LlmsEntity(
        LLMType.Azure,
        LLMType.Azure,
        'Wrapper around OpenAI Chat Large language models.',
        'Wrapper around OpenAI Chat Large language models.',
        [],
        [
          {
            label: 'api_key',
            key: LLMFields.API_KEY,
            hint: 'The API key for the Azure API.',
            type: NodeFieldType.TextField,
            required: true,
            value: process.env.OPENAI_API_KEY || '', // default value
            accept: {},
          },
          {
            label: 'api_base',
            key: LLMFields.API_BASE,
            hint: 'The base URL for the Azure API.',
            type: NodeFieldType.TextField,
            required: true,
            value: process.env.OPENAI_API_BASE || '', // default value
            accept: {},
          },
          {
            label: 'api_engine',
            key: LLMFields.API_ENGINE,
            hint: 'The engine for the Azure API.',
            type: NodeFieldType.TextField,
            required: true,
            value: process.env.OPENAI_API_ENGINE || '', // default value
            accept: {},
          },
        ],
        [],
        {
          type: NodeTypes.LLMs,
          name: NodeTypes.LLMs,
          uuid: '',
          functionName: '',
          functionDesc: '',
          required: false,
        },
        '0.0.1',
      ),
    ];
  }
  async getDeviceTemplates() {
    return [
      new DevicesEntity(
        DeviceType.OT2,
        DeviceType.OT2,
        'Description:\n' +
          '\n' +
          '\n' +
          '●OT-2 is a laboratory automation robot designed for high-precision liquid handling and sample analysis.\n' +
          '●It is capable of autonomously performing a variety of lab tasks such as pipetting, mixing, centrifugation, and temperature control.\n' +
          '●Supports a range of experimental container formats, including test tubes, petri dishes, and microplates.\n' +
          '\n' +
          'Interface and Communication:\n' +
          '\n' +
          '\n' +
          '●OT-2 communicates with external systems using the TCP/IP protocol.\n' +
          '●Supports RESTful API, allowing commands to be sent and responses received via HTTP requests.\n' +
          '●Provides command and response formats in JSON, ensuring cross-platform compatibility.\n',
        [],
        [
          {
            label: 'api_base',
            key: DeviceFields.API_BASE,
            hint: 'The base URL for the OT2 API.',
            type: NodeFieldType.TextField,
            required: true,
            value: 'http://127.0.0.1:31950', // default value
            accept: {},
          },
          {
            label: 'jetson_api_base',
            key: DeviceFields.JETSON_API_BASE,
            hint: 'The base URL for the Jetson API.',
            type: NodeFieldType.TextField,
            required: true,
            value: 'http://127.0.0.1:8081', // default value
            accept: {},
          }, // Jetson
        ],
        [],
        {
          type: NodeTypes.Devices,
          name: NodeTypes.Devices,
          uuid: '',
          functionName: '',
          functionDesc: '',
          required: false,
        },
        '0.0.1',
      ),
    ];
  }

  async getPromptTemplates() {
    return [
      new PromptEntity(
        PromptType.PromptTemplate,
        PromptType.PromptTemplate,
        'Schema to represent a prompt for an LLM',
        'Schema to represent a prompt for an LLM',
        [],
        [
          {
            label: 'Prompt',
            key: PromptFields.PROMPT,
            hint: 'The prompt for the LLM.',
            required: true,
            accept: {},
            type: NodeFieldType.TextArea,
          },
          {
            label: 'Role',
            key: PromptFields.ROLE,
            hint: 'The role for the LLM.',
            required: true,
            accept: {},
            type: NodeFieldType.TextField,
          },
        ],
        [],
        {
          type: NodeTypes.Prompts,
          name: NodeTypes.Prompts,
          uuid: '',
          functionName: '',
          functionDesc: '',
          required: false,
        },
        '0.0.1',
      ),
    ];
  }

  async getInputTemplates() {
    return [
      new InputEntity(
        InputType.InputNode,
        InputType.InputNode,
        'Input node, used for automatic docking of inputs',
        [],
        [
          {
            label: 'Input',
            key: InputFields.INPUT,
            hint: 'Input...',
            value: '',
            required: true,
            accept: {},
            type: NodeFieldType.Inputs,
          },
        ],
        [],
        {
          type: NodeTypes.Input,
          name: NodeTypes.Input,
          uuid: '',
          functionName: '',
          functionDesc: '',
          required: false,
        },
        '0.0.1',
      ),
    ];
  }

  async getToolTemplates() {
    return [
      new ToolEntity(
        ToolType.GOOGLE_SEARCH,
        ToolType.GOOGLE_SEARCH,
        'A low-cost Google Search APl.Useful for when youneed to answer questions about biomedical health-related.Input should be a search query.',
        [],
        [
          {
            label: 'Google Search Api Key',
            key: ToolFields.GOOGLE_API_KEY,
            hint: 'google api key...',
            value: process.env.GOOGLE_SEARCH_API_KEY || '',
            required: true,
            accept: {},
            type: NodeFieldType.Inputs,
          },
          {
            label: 'Google Search Engine ID',
            key: ToolFields.GOOGLE_ENGINE_ID,
            hint: 'google engine id...',
            value: process.env.GOOGLE_SEARCH_API_CX || '',
            required: true,
            accept: {},
            type: NodeFieldType.Inputs,
          },
        ],
        [],
        {
          type: NodeTypes.Tools,
          name: NodeTypes.Tools,
          uuid: '',
          functionName: '',
          functionDesc: '',
          required: false,
        },
        '0.0.1',
      ),

      new ToolEntity(
        ToolType.WIKIPEDIA_SEARCH,
        ToolType.WIKIPEDIA_SEARCH,
        'A tool for interacting with and fetching data from the Wikipedia API.',
        [],
        [
          {
            label: 'Wikipedia Top K Results',
            key: ToolFields.WIKIPEDIA_TOP_K_RESULTS,
            hint: 'wikipedia top k results...',
            value: '3',
            required: true,
            accept: {},
            type: NodeFieldType.Inputs,
          },
          {
            label: 'Wikipedia Max Doc Content Length',
            key: ToolFields.WIKIPEDIA_MAX_DOC_CONTENT_LENGTH,
            hint: 'wikipedia max doc content length...',
            value: '4000',
            required: true,
            accept: {},
            type: NodeFieldType.Inputs,
          },
        ],
        [],
        {
          type: NodeTypes.Tools,
          name: NodeTypes.Tools,
          uuid: '',
          functionName: '',
          functionDesc: '',
          required: false,
        },
        '0.0.1',
      ),
    ];
  }
}
