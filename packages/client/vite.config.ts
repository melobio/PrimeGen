import { fileURLToPath, URL } from 'node:url'

import { defineConfig } from 'vite'
import vue from '@vitejs/plugin-vue'
import vueJsx from '@vitejs/plugin-vue-jsx'
import vuetify from 'vite-plugin-vuetify'

// https://vitejs.dev/config/
export default defineConfig({
  base: '/xpcr',
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
      "/xpcr-api": {
        target: "http://localhost:3001",
        changeOrigin: true,
        rewrite: (path: string) => {
          // console.log(path)
          return path
        },
      },
      "/xpcr-ws": {
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
      "/xpcr-search-api": {
        target: "http://localhost:8081",
        changeOrigin: true,
        rewrite: (path) =>{
          const newPath =  path.replace(/^\/xpcr-search-api/, '');
          return newPath
        }
      },
      "/xpcr-search-dev": {
        target: "https://chat.mgi-tech.com/",
        changeOrigin: true,
        rewrite: (path) =>{
          return path
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
