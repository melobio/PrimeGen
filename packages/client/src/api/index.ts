import axios from 'axios'
import Cookies from 'js-cookie'
import { XResponse, type AdapterRes } from '@/api/response'
import type { ToolsName } from '@/pages/main/tools/ToolItem'
import type { AgentNode } from "@xpcr/common/src/agent";
import type {DeviceNode, InputNode, LLMNode, PromptNode, Node, ToolNode} from "@xpcr/common";

const client = axios.create({
  timeout: 10000,
  baseURL: import.meta.env.VITE_APP_API
})

const client_adpater = axios.create({
  timeout: 60000,
  baseURL: import.meta.env.VITE_APP_ADAPTER_API
})

const client_search = axios.create({
  timeout: 60000,
  baseURL: import.meta.env.VITE_APP_SEARCH_API
})

client.interceptors.request.use((request) => {
  if (Cookies.get('x-token')) {
    request.headers['Authorization'] = 'Bearer ' + Cookies.get('x-token')
  }
  return request
})
client.interceptors.response.use((response) => {
  return response.data;
}, () => {
  return new XResponse(0, '未知错误', false, null);
})

client_adpater.interceptors.response.use((response): any => {

  const adapterRes: AdapterRes = response.data;
  if (adapterRes.code === 200) {
    return new XResponse(adapterRes.code, adapterRes.msg, true, adapterRes.data);
  }
  return new XResponse(adapterRes.code, adapterRes.msg, false, null);
}, () => {
  return new XResponse(0, '未知错误', false, null);
})

export function hello() {
  return client.get<string>('/hello')
}
export function login<T>(
  data: {
    loginType: 'userName'| 'email' | 'phone';
    userName?: string;
    password?: string;
    phone?: string;
    email?: string;
    smsCode?: string;
  },
): Promise<XResponse<T>>{
  return client.post('/auth/login', data);
}

export function experiments<T>(
): Promise<XResponse<T>> {
  return client.get('/experiments');
}

export function addConWithExp<T>(
  ): Promise<XResponse<T>> {
    return client.post('/experiments/addConWithExp');
  }

export function chats<T>(
): Promise<XResponse<T>> {
  return client.get('/conversations');
}

export function messages<T>(
  conversationUUID: string,
): Promise<XResponse<T>> {
  return client.get(`conversations/${conversationUUID}/messages`);
}

export function deleteMessageById<T>(
  messageEntity: any,
): Promise<XResponse<T>> {
  return client.post(`conversations/deleteMessageById`, messageEntity);
}

export function deleteConversationByUUID<T>(
  conversationUUID: string,
): Promise<XResponse<T>> {
  return client.post(`conversations/${conversationUUID}/del`);
}

export function toolsAgents<T>(
): Promise<XResponse<T>> {
  return client.get('/tools/agents');
}
export function toolsAgentTemplates<T>(
): Promise<XResponse<T>> {
  return client.get('/tools/template/agents');
}
export function updateAgents(
  uuid: string,
  data: AgentNode,
): Promise<XResponse<any>> {
  return client.post(`/tools/agents/${uuid}`, data);
}
export function toolsDevices<T>(
): Promise<XResponse<T>> {
  return client.get('/tools/devices');
}
export function toolsDeviceTemplates<T>(
): Promise<XResponse<T>> {
  return client.get('/tools/template/devices');
}
export function updateDevice(
  uuid: string,
  data: DeviceNode,
): Promise<XResponse<any>> {
  return client.post(`/tools/devices/${uuid}`, data);
}
export function checkDeviceState(uuid: string): Promise<XResponse<any>> {
  return client.get(`/tools/devices/${uuid}/state`);
}
export function toolsInputs<T>(
): Promise<XResponse<T>> {
  return client.get('/tools/inputs');
}
export function toolsInputTemplates<T>(
): Promise<XResponse<T>> {
  return client.get('/tools/template/inputs');
}

export function updateInput(
  uuid: string,
  data: InputNode,
): Promise<XResponse<any>> {
  return client.post(`/tools/inputs/${uuid}`, data);
}
export function toolsPrompts<T>(
): Promise<XResponse<T>> {
  return client.get('/tools/prompts');
}
export function toolsPromptTemplates<T>(
): Promise<XResponse<T>> {
  return client.get('/tools/template/prompts');
}
export function updatePrompt(
  uuid: string,
  data: PromptNode
): Promise<XResponse<any>> {
  return client.post(`/tools/prompts/${uuid}`, data);
}

export function toolsTools<T>(
): Promise<XResponse<T>> {
  return client.get('/tools/tools');
}
export function toolsToolTemplates<T>(
): Promise<XResponse<T>> {
  return client.get('/tools/template/tools');
}
export function updateTool(
  uuid: string,
  data: ToolNode
): Promise<XResponse<any>> {
  return client.post(`/tools/tools/${uuid}`, data);
}

export function toolsLLms<T>(
): Promise<XResponse<T>> {
  return client.get('/tools/llms');
}
export function toolsLLMTemplates<T>(
): Promise<XResponse<T>> {
  return client.get('/tools/template/llms');
}
export function updateLLMs(
  uuid: string,
  data: LLMNode,
): Promise<XResponse<any>> {
  return client.post(`/tools/llms/${uuid}`, data);
}

export function updateLLms(
  uuid: string,
  data: LLMNode
): Promise<XResponse<any>> {
  return client.post(`/tools/llms/${uuid}`, data);
}

export function updateCommonToolItem(
  uuid: string,
  type: ToolsName,
  data: Node,
): Promise<XResponse<any>> {
  return client.post(`/tools/${type.toLowerCase()}/${uuid}`, data);
}

export function commonToolInfo<T>(
  uuid: string,
  type: ToolsName,
): Promise<XResponse<T>> {
  return client.get(`/tools/${type.toLowerCase()}/${uuid}`);
}

export function uploadPrimerExcel<T>(formData: FormData): Promise<XResponse<T>> {
  return client_adpater.post(`/file/upload`, formData, {
    headers: {
      'Content-Type': 'multipart/form-data',
    },
  });
}

export function uploadSequenceFiles(formData: FormData){
  return client_search.post(`/upload`, formData, {
    headers: {
      'Content-Type': 'multipart/form-data',
    },
  });
}