import { fileURLToPath, URL } from 'node:url'

import { defineConfig } from 'vite'
import vue from '@vitejs/plugin-vue'
import vueJsx from '@vitejs/plugin-vue-jsx'
import vuetify from 'vite-plugin-vuetify'

// https://vitejs.dev/config/
export default defineConfig({
  plugins: [
    vue(),
    vueJsx(),
    vuetify({ autoImport: true })
  ],
  resolve: {
    alias: {
      '@': fileURLToPath(new URL('./src', import.meta.url))
    }
  },
  server: {
    proxy: {
      "/api": {
        target: "http://localhost:3001",
        changeOrigin: true,
        rewrite: (path: string) => {
          // console.log(path)
          return path
        },
      },
      "/pcr-ws": {
        target: "ws://localhost:3001",
        changeOrigin: true,
        ws: true,
        rewrite: (path: string) => {
          // console.log(path)
          return path
        },
      },
      "/v1/files": {
        target: "http://localhost:3001",
        changeOrigin: true,
        rewrite: (path: string) => {
          // console.log(path)
          return path
        },
      },
      "/primer": {
        target: "http://172.16.225.244:9516",
        changeOrigin: true,
        rewrite: (path: string) => {
          // console.log('path=', path)
          return path
        }
      },
      "/search": {
        target: "http://172.16.47.11:8081",
        changeOrigin: true,
        rewrite: (path) =>{
          const newPath =  path.replace(/^\/search/, '');
          return newPath
        }
      }
    }
  },
  build: {
    rollupOptions: {
      external: ['build/labware.json']
    }
  }
})
