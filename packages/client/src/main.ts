import './assets/main.css'

let process: any;

import { createApp } from 'vue'
import { createPinia } from 'pinia'
import vuetify from '@/plugins/vuetify'
import capitalize from 'lodash/capitalize'
import startCase from 'lodash/startCase'
import { titleCase } from '@xpcr/shared/src/titleCase'
import i18next from "i18next";
import I18NextVue from "i18next-vue";
import { resources } from '@/assets/localization'
i18next.init({
  resources,
  lng: 'en',
  fallbackLng: 'en',
  debug: process.env.NODE_ENV === 'development',
  ns: [
    'shared',
    'robot_advanced_settings',
    'robot_calibration',
    'robot_connection',
    'robot_controls',
    'robot_info',
    'top_navigation',
  ],
  defaultNS: 'shared',
  interpolation: {
    escapeValue: false, // not needed for react as it escapes by default
    format: function (value, format, lng) {
      if (format === 'upperCase') return value.toUpperCase()
      if (format === 'capitalize') return capitalize(value)
      if (format === 'sentenceCase') return startCase(value)
      if (format === 'titleCase') return titleCase(value)
      return value
    },
  },
  keySeparator: false, // use namespaces and context instead
  saveMissing: true,
  missingKeyHandler: (lng, ns, key) => {
    process.env.NODE_ENV === 'test'
      ? console.error(`Missing ${lng} Translation: key={${key}} ns={${ns}}`)
      : console.warn(`Missing ${lng} Translation: key={${key}} ns={${ns}}`)
  },
})

import App from './App.vue'
import router from './router'

const app = createApp(App)

app.use(createPinia())
app.use(vuetify)
app.use(router)
app.use(I18NextVue, { i18next })
// 去掉警告消息
app.config.warnHandler = () => null;

app.mount('#app')
