<template>
  <div class='local'>
    <v-avatar :image='User' size='48'></v-avatar>
    <div class='arrow'>
    </div>
    <div v-html="markedContent" class="content">
    </div>
  </div>
</template>

<script setup lang='ts'>
import '@/styles/markdown.css'
import 'highlight.js/styles/atom-one-dark-reasonable.min.css';
import { computed } from "vue";
import User from '@/assets/user.png'
import mdKatex from '@traptitech/markdown-it-katex'
import MarkdownIt from 'markdown-it'
import mila from 'markdown-it-link-attributes'
import hljs from 'highlight.js';
import type {Message} from "@/stores/chat-types";
const props = defineProps<Message>()

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
const highlightBlock = (str: string, lang?: string) => {
  return `<pre class="code-block-wrapper"><div class="code-block-header"><span class="code-block-header__lang">${lang}</span>${lang && '|'}<span class="code-block-header__copy">copy</span></div><code class="hljs code-block-body ${lang}">${str}</code></pre>`
}
const capitalizeFirstLetterRegex = (str:string) => {
    return str.replace(/^\w/, (match) => match.toUpperCase());
}
mdi.use(mila, { attrs: { target: '_blank', rel: 'noopener' } })
mdi.use(mdKatex, { blockClass: 'katexmath-block', errorColor: ' #cc0000' })
  const markedContent = computed(() => {
  let htmlContent = mdi.render(props.content)
  return htmlContent;
})
</script>

<style scoped lang='scss'>
.local {
  display: flex;
  flex-direction: row-reverse;
  align-items: end;
  margin-top: 25px;
  padding-left: 83px;
  padding-right: 20px;
  .arrow {
    width: 15px;
    height: 15px;
    flex-shrink: 0;
    background-color: #31DA9F;
    opacity: 0.8;
    &::after {
      width: 15px;
      height: 15px;
      content: '';
      background-color: #1E2032;
      float: left;
      border-bottom-left-radius: 15px;
    }
  }
  .content {
    display: flex;
    align-items: center;
    background-color: #31DA9F;
    opacity: 0.8;
    border-top-right-radius: 10px;
    border-top-left-radius: 10px;
    border-bottom-left-radius: 10px;
    color: #101010;
    max-width: 600px;
    min-width: 40px;
    padding: 8px 10px;
    min-height: 40px;
    white-space: pre-wrap;
  
    }
}
</style>