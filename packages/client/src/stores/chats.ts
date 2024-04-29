import { defineStore } from 'pinia'
import { ref, watch } from 'vue'
import type { Chat, Message, ConversationUUIDtoPrimerFileID, GeneticDisorderInfo, CancerOptionInfo, ProteinMutationInfo, PathogenDrugInfo,SpeciesIdentificationInfo } from '@/stores/chat-types'
import { chats as getChats, messages as getMessages , deleteConversationByUUID} from '@/api';
import { io, Socket } from "socket.io-client";
import * as _ from 'lodash';
import {
  TEXT_ANSWER_CREATE,
  TEXT_ANSWER_DONE,
  TEXT_ANSWER_GENERATING,
  TEXT_QUESTION,
  JETSON_BAD_TIP,
  UPLOAD_PRIMER_EXCEL,
  VOICE_QUESTION,
  SEND_OPTION_NCBI_SEARCH,
} from "@/stores/chat-types";

export interface BadTipInfo {
  runId: string,
  bad_tip: number,
  pcr: number,
  tip: number,
  conversationUUID: string,
  data: string,
  plot: Array<number>
}

interface GeneratingMessage {
  success: boolean,
  data: Message,
  finishReason: string,
  isCreate: boolean,
  conversationUUID: string,
}

function initSocket() {
  // const socketOptions = {
  //   path: '/pcr-ws',
  //   autoConnect: false,
  //   transportOptions: {
  //     polling: {
  //       extraHeaders: {
  //         Authorization: Cookies.get('x-token'),
  //       }
  //     }
  //   }
  // }
  // console.log('initSocket', Cookies.get('x-token'));
  return io('/', {
    path: import.meta.env.VITE_APP_WS_PATH,
    autoConnect: true,
  });
}

export const useChatStore = defineStore('chat', () => {
  const socket = ref<Socket>();
  const chats = ref<Chat[]>([]);
  const currentChat = ref<Chat>({} as Chat);
  const currentChatMessages = ref<Message[]>([]);
  const isRecordingAudio = ref(false);
  const badTipInfo = ref<BadTipInfo>({
    runId: '',
    bad_tip: 0,
    pcr: 0,
    tip: 0,
    conversationUUID: '',
    data: '',
    plot: []
  });

  const generatingMessages: GeneratingMessage[] = [];
  // const isGenerating = ref(false);
  const cvstoexper = ref<ConversationUUIDtoPrimerFileID[]>([])

  // 消息更新时，是否滚动到底部
  const scrollBottom = ref(false);
  const shouldScrollToBottomWhileGenerating = ref(false);

  watch(() => currentChat.value, (newVal) => {
    if (newVal?.uuid) {
      getAllMessages(newVal.uuid).then();
    }
  });

  async function init() {
    socket.value?.disconnect();
    socket.value = initSocket();
    socket.value.on(TEXT_ANSWER_GENERATING, async (generatingMessage: GeneratingMessage) => {
      if (generatingMessage.success) {
        generatingMessages.push(generatingMessage);
        if (generatingMessage.finishReason == 'stop') {
          scrollBottom.value = true;
          await getAllMessages(currentChat.value.uuid);        
        }else{
          updateGenerating(generatingMessage.data.uuid,true,false)
        }
      }
    });
    socket.value.on(TEXT_ANSWER_CREATE, (generatingMessage: GeneratingMessage) => {
      if (generatingMessage.success) {
        generatingMessages.push({ ...generatingMessage, isCreate: true });
      }
    });
    socket.value.on(TEXT_ANSWER_DONE,async ({ conversationUUID }: { conversationUUID: string }) => {
      console.log('TEXT_ANSWER_DONE', conversationUUID);
      await getAllChats();
      await getAllMessages(currentChat.value.uuid);        
      updateGenerating(currentChat.value.uuid,false,true)
    })
    socket.value.on(JETSON_BAD_TIP, ({ conversationUUID, bad_tip, pcr, tip, data, plot, runId }: BadTipInfo) => {
      setBadTipInfo({ runId, conversationUUID, bad_tip, pcr, tip, data, plot });
    })

    socket.value.on(UPLOAD_PRIMER_EXCEL, (generatingMessage: GeneratingMessage) => {
      console.log('recieve message, ready to upload primer excel');
      console.log(generatingMessage)
      currentChatMessages.value.push({ ...generatingMessage.data, generating: false });
    })
    startRefreshGeneratingMessages();
  }
  function setRecordingAudioState(val: boolean) {
    isRecordingAudio.value = val;
  }
  function setBadTipInfo(info: BadTipInfo) {
    badTipInfo.value = info;
  }
  // 每50ms检查一次，是否有新的生成消息，有的话就显示出来
  // 目的是为了让用户看到生成的过程
  function startRefreshGeneratingMessages() {
    setInterval(() => {
      if (generatingMessages.length) {
        const generatingMessage = generatingMessages.shift();
        if (generatingMessage) {
          const targetMessage = currentChatMessages.value.find((item) => item.id === generatingMessage.data.id);
          if (targetMessage) {
            if (generatingMessage.finishReason === 'stop') {
              _.assign(targetMessage, { generating: false });
            } else {
              _.assign(targetMessage, {
                content: generatingMessage.data.content,
                header: generatingMessage.data.header,
                protocolAnalysis: generatingMessage.data.protocolAnalysis,
                currentCommandIndex: generatingMessage.data.currentCommandIndex,
                faultCheckResult: generatingMessage.data.faultCheckResult,
                agentType: generatingMessage.data.agentType,
                optionInfo: generatingMessage.data.optionInfo,
              });
            }
          } else if (generatingMessage.isCreate && generatingMessage.data.conversationUUID == currentChat.value.uuid) {
            currentChatMessages.value.push({ ...generatingMessage.data, generating: true });
          }
          if (shouldScrollToBottomWhileGenerating.value) {
            scrollBottom.value = true;
          }
        }
      }
    }, 50);
  }

 async function updateGenerating(uuid: string,isGenerating:boolean,setCurrent:boolean){
  if(setCurrent){currentChat.value.isGenerating = isGenerating;}
    chats.value.forEach((item:Chat,index:number)=>{
      if(item.uuid == uuid){
        chats.value[index].isGenerating = isGenerating;
      }
    });
  }
  async function getAllChats() {
    const { data = [], success } = await getChats<Chat[]>()
    if (success) {
      const newChats = data.map(newItem=>{
        let newChatItem:Chat |null = null;
        chats.value.forEach( subItem => {
          if(subItem.id == newItem.id){
            newChatItem = {...subItem,...newItem}
          }
        })
        if(!newChatItem){
          newChatItem = { ...newItem,isGenerating:false,isRecordingAudio:false }
        }
        return newChatItem;
      })
      chats.value = newChats;
      if (data?.length && !currentChat.value?.uuid) {
        currentChat.value = data[0];
      }
    }
  }
  async function getAllMessages(conversationUUID: string) {
    getMessages<Message[]>(conversationUUID).then(({ data, success }) => {
      if (success) {
        data?.forEach((item) => item.generating = false);
        currentChatMessages.value = data || [];
        scrollBottom.value = true;
      }
    });
  }
  // AI
  async function sendText(text: string) {
    if (!currentChat.value?.uuid) {
      console.warn('currentChat empty');
      return;
    }
    if (!socket.value || socket.value?.connected === false) {
      console.warn('socket empty or disconnected');
      return;
    }
    updateGenerating(currentChat.value.uuid,true,true)
    const { success, data }: { success: boolean, data: Message[] } =
      await socket.value.emitWithAck(TEXT_QUESTION, { text, conversationUUID: currentChat.value?.uuid });

    if (success) {
      data[data.length - 1].generating = false; // 正在生成中
      currentChatMessages.value = [...currentChatMessages.value, ...data];
      scrollBottom.value = true;
    }
  }

  async function sendVoice(voice: Blob) {
    if (!currentChat.value?.uuid) {
      console.warn('currentChat empty');
      return;
    }
    if (!socket.value || socket.value?.connected === false) {
      console.warn('socket empty or disconnected');
      return;
    }
    // blob to base64
    const base64Voice = await blobToBase64(voice);

    console.log('base64Voice', base64Voice.substring(0, 100) + '...');
    const [_, base64] = base64Voice.split(';base64,');
    updateGenerating(currentChat.value.uuid,true,true)
    const { success, data }: { success: boolean, data: Message[] } =
      await socket.value.emitWithAck(VOICE_QUESTION, { voice: base64, conversationUUID: currentChat.value?.uuid });

    if (success) {
      data[data.length - 2].generating = false; // 用户发出的直接标记为生成结束
      data[data.length - 1].generating = true; // 正在生成中
      currentChatMessages.value = [...currentChatMessages.value, ...data];
      scrollBottom.value = true;
    }
  }
  async function handleDeleteConversationByUUID(uuid: string){
    const res =  await deleteConversationByUUID(uuid)
    if(res.code == 200){
      await getAllChats();
      if(currentChat.value.uuid === uuid){
        currentChat.value = chats.value[chats.value.length-1]
      }
    }
  }

  async function sendSelectOptionToNcbiSearch(
    text: string,
    info: ProteinMutationInfo| CancerOptionInfo | GeneticDisorderInfo | PathogenDrugInfo | SpeciesIdentificationInfo,
    message: Message) {
    if (!currentChat.value?.uuid) {
      console.warn('currentChat empty');
      return;
    }
    if (!socket.value || socket.value?.connected === false) {
      console.warn('socket empty or disconnected');
      return;
    }
    updateGenerating(currentChat.value.uuid,true,true)
    let res:{ success: boolean, data: Message[] }= {success:false,data:[]};
      res= await socket.value.emitWithAck(SEND_OPTION_NCBI_SEARCH, {
          text,
          conversationUUID: currentChat.value?.uuid,
          option: info,
          messageId: message.id,
        });
    if (res.success) {
      currentChatMessages.value = [...currentChatMessages.value, ...res.data];
      scrollBottom.value = true;
    }
  }

  async function blobToBase64(blob: Blob) {
    return new Promise<string>((resolve) => {
      const reader = new FileReader();
      reader.readAsDataURL(blob);
      reader.onload = () => {
        const base64 = reader.result?.toString();
        if (base64) {
          resolve(base64);
        }
      }
    });
  }

  async function sendOption(option: 'stop' | 'play', runId: string) {
    if (!currentChat.value?.uuid) {
      console.warn('currentChat empty');
      return;
    }
    if (!socket.value || socket.value?.connected === false) {
      console.warn('socket empty or disconnected');
      return;
    }
    const { success, data }: { success: boolean, data: '' } =
      await socket.value.emitWithAck(JETSON_BAD_TIP, { runId, option, conversationUUID: currentChat.value?.uuid });
    return { success, data };
  }

  function logout() {
    chats.value = [];
    currentChat.value = {} as Chat;
    currentChatMessages.value = [];
    generatingMessages.length = 0;
  }

  return {
    init,
    socket,
    chats,
    currentChat,
    currentChatMessages,
    isRecordingAudio,
    badTipInfo,
    getAllChats,
    getAllMessages,
    sendText,
    sendOption,
    handleDeleteConversationByUUID,
    sendSelectOptionToNcbiSearch,
    sendVoice,
    // updateUploadSequenceFileMessage,
    scrollBottom,
    shouldScrollToBottomWhileGenerating,
    logout,
    setRecordingAudioState,
    setBadTipInfo,
    // isGenerating,
    cvstoexper
  }
});