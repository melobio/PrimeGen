<template>
  <div class='remote'>
    <v-avatar :image='remoteIcon' size='48'></v-avatar>
    <div class='arrow'>
    </div>
    <div class='content'>
      <!-- <span v-html="header" v-if="!uploadFlag"></span> -->
      <JetFaultCheckResult v-for="(item, index) in globalFailCheckResult" :key="index" v-bind="item">
      </JetFaultCheckResult>
      <OTCommands v-if="protocolAnalysis" :current-command-index="currentCommandIndex == null ? -1 : currentCommandIndex"
        :command-fail-check-result="commandFailCheckResult" :protocol-analysis="protocolAnalysis" />
      <div class="markdown-body" ref="textRef">
        <div v-html="markedContent" v-if="markedContent && !uploadFlag" class="pre text-content"
          :class="{ 'markdown-body-generate': (generating && markedContent) }"></div>
        <!-- 思考中的光标 -->
        <div v-if="generating && (!markedContent)" class="text-cursor" />
        <!-- 物种鉴定 -->
        <SpeciesIdentification v-if="showSpeciesIdentification" v-bind="props"/>
        <!--遗传病基因选项 -->
        <GeneticDiseasesOptions v-if="showGeneticDiseasesOptions" v-bind="props"/>
        <!--cancer基因选项-->
        <CancerInfoOptions v-if="showCancerInfoOptions"  v-bind="props"/>
        <!--蛋白质突变选项 或者 用户直接下载基因序列 -->
        <ProteinMutationOptions v-if="showProteinMutationOptions" v-bind="props"/>
        <!--病原体药物抗性选项 -->
        <PathogenDrugOptions v-if="showPathogenDrugOptions" v-bind="props"/>
        <UploadPrimer v-if="uploadFlag" v-bind="props"/>
        <InitiativeStartPrimerDesign v-if="showInitiativeStartPrimerDesign" v-bind="props"/>
      </div>
    </div>
  </div>
</template>

<script setup lang='ts'>
import '@/styles/markdown.css'
import MarkdownIt from 'markdown-it'
import mdKatex from '@traptitech/markdown-it-katex'
import hljs from 'highlight.js';
import 'highlight.js/styles/atom-one-dark-reasonable.min.css';
import mila from 'markdown-it-link-attributes'
import GeneticDiseasesOptions from './GeneticDiseasesOptions.vue'
import { copyToClip } from '@/utils/copy'
import IconPlanner from '@/assets/planner.png';
import IconFaultAgent from '@/assets/fault_agent.png';
import IconInternetSearchAgent from '@/assets/internet_search_agent.png';
import IconSequenceSearchAgent from '@/assets/sequence_search_agent.png';
import IconCodeExecutionAgent from '@/assets/code_execution_agent.png';
import IconPrimerDesignAgent from '@/assets/primer_design_agent.png';
import { UPLOAD_PRIMER_EXCEL, NCBIFunctionType,INITIATIVE_START_PRIMER_DESIGN,type Message } from "@/stores/chat-types";
import { computed, onMounted, onUnmounted, onUpdated, reactive, ref, watchEffect,watch } from "vue";
import OTCommands from "@/pages/main/chat/chat-detail/OTCommands.vue";
import JetFaultCheckResult from "@/pages/main/chat/chat-detail/JetFaultCheckResult.vue";
import { AgentType } from "@xpcr/common/src/agent";
import UploadPrimer from '@/pages/main/chat/chat-detail/UploadPrimer.vue'
import CancerInfoOptions from './CancerInfoOptions.vue'
import ProteinMutationOptions from './ProteinMutationOptions.vue'
import PathogenDrugOptions from './PathogenDrugOptions.vue'
import SpeciesIdentification from './SpeciesIdentification.vue'
import InitiativeStartPrimerDesign from './InitiativeStartPrimerDesign.vue'
const props = defineProps<Message>()
const textRef = ref<HTMLElement>()

let uploadFlag = ref<Boolean>(false)

watchEffect(() => {
  if (props.header === UPLOAD_PRIMER_EXCEL) {
    uploadFlag.value = true;
  } else {
    uploadFlag.value = false;
  }
});
const showInitiativeStartPrimerDesign = computed(()=>{
  return props.header == INITIATIVE_START_PRIMER_DESIGN
})
const showProteinMutationOptions = computed(()=>{
  let searchType = props?.optionInfo?.search_type;
  return (searchType == NCBIFunctionType.protein_mutation_type||searchType == NCBIFunctionType.download_type) && props.agentType == AgentType.SEQUENCE_SEARCH ;
})

const showPathogenDrugOptions = computed(()=>{
  return props?.optionInfo?.search_type == NCBIFunctionType.pathogen_drug_resistance_type && props.agentType == AgentType.SEQUENCE_SEARCH ;
})

const showGeneticDiseasesOptions = computed(()=>{
  return props?.optionInfo?.search_type == NCBIFunctionType.genetic_disorder_type && props.agentType == AgentType.SEQUENCE_SEARCH && props?.optionInfo?.options.length>0;
})
const showCancerInfoOptions = computed(()=>{
  return props?.optionInfo?.search_type == NCBIFunctionType.cancer_type && props.agentType == AgentType.SEQUENCE_SEARCH;
})

const showSpeciesIdentification = computed(()=>{
  return props?.optionInfo?.search_type == NCBIFunctionType.species_identification_type;
})

const globalFailCheckResult = computed(() => {
  return props.faultCheckResult?.filter(item => !item.commandIndex) ?? [];
});

const commandFailCheckResult = computed(() => {
  return props.faultCheckResult?.filter(item => item.commandIndex) ?? [];
});
const remoteIcon = computed(() => {
  switch (props.agentType) {
    case AgentType.FAULT:
      return IconFaultAgent;
    case AgentType.INTERNET_SEARCH:
      return IconInternetSearchAgent;
    case AgentType.SEQUENCE_SEARCH:
      return IconSequenceSearchAgent;
    case AgentType.CODE_EXECUTION:
      return IconCodeExecutionAgent;
    case AgentType.PRIMER_DESIGN:
      return IconPrimerDesignAgent;
    default:
      return IconPlanner;
  }
});
const markedContent = computed(() => {
  let htmlContent = mdi.render(props.content)
  return htmlContent;
})

const addCopyEvents = () => {
  if (textRef.value) {
    const copyBtn = textRef.value.querySelectorAll('.code-block-header__copy')
    copyBtn.forEach((btn) => {
      btn.addEventListener('click', () => {
        const code = btn.parentElement?.nextElementSibling?.textContent
        if (code) {
          copyToClip(code).then(() => {
            btn.textContent = `已复制!`
            setTimeout(() => {
              btn.textContent = `复制`
            }, 1000)
          })
        }
      })
    })
  }
}

const removeCopyEvents = () => {
  if (textRef.value) {
    const copyBtn = textRef.value.querySelectorAll('.code-block-header__copy')
    copyBtn.forEach((btn) => {
      btn.removeEventListener('click', () => { })
    })
  }
}

onMounted(() => {
  addCopyEvents()
})
onUpdated(() => {
  addCopyEvents()
})
onUnmounted(() => {
  removeCopyEvents()
})
const mdi = new MarkdownIt({
  html: false,
  linkify: true,
  highlight(code, language) {
    const validLang = !!(language && hljs.getLanguage(language))
    if (validLang) {
      const lang = capitalizeFirstLetterRegex(language) ?? ''
      return highlightBlock(hljs.highlight(code, { language: lang }).value, lang)
    }
    return highlightBlock(hljs.highlightAuto(code).value, '')
  },
})

mdi.use(mila, { attrs: { target: '_blank', rel: 'noopener' } })
mdi.use(mdKatex, { blockClass: 'katexmath-block', errorColor: ' #cc0000' })
const highlightBlock = (str: string, lang?: string) => {
  return `<pre class="code-block-wrapper"><div class="code-block-header"><span class="code-block-header__lang">${lang}</span>${lang && '|'}<span class="code-block-header__copy">复制</span></div><code class="hljs code-block-body ${lang}">${str}</code></pre>`
}
const capitalizeFirstLetterRegex = (str:string) => {
    return str.replace(/^\w/, (match) => match.toUpperCase());
}
</script>

<style  lang='scss'>
.remote {
  display: flex;
  align-items: end;
  margin-top: 25px;
  padding-left: 20px;
  padding-right: 83px;

  .arrow {
    width: 15px;
    height: 15px;
    flex-shrink: 0;
    background-color: rgb(52, 54, 72);

    &::after {
      width: 15px;
      height: 15px;
      content: '';
      background-color: #1E2032;
      float: left;
      border-bottom-right-radius: 15px;
    }
  }

  .content {
    background-color: rgb(52, 54, 72);
    border-top-right-radius: 10px;
    border-top-left-radius: 10px;
    border-bottom-right-radius: 10px;
    color: #ccc;
    max-width: 630px;
    min-width: 40px;
    padding: 8px 10px;
    min-height: 40px;
    white-space: pre-wrap;

    .text-content {
      word-break: break-all;
      list-style-position: inside;
      display: inline-flex;
      flex-direction: column;
      margin-top: 0;
      margin-bottom: 0;
      font-size: 15px;
      word-wrap: normal;
    }

    .text-cursor:after {
      content: '';
      font-weight: 700;
      width: 2px;
      height: 18px;
      margin-left: 6px;
      vertical-align: baseline;
      background-color: #1484FC;
      display: inline-block;
      animation: blink 1s linear infinite;
    }

    @keyframes blink {

      0%,
      100% {
        opacity: 1;
      }

      50% {
        opacity: 0;
      }
    }
  }
}</style>import type { fromPairs } from 'lodash';
