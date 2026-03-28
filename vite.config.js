import { defineConfig } from 'vite'
import react from '@vitejs/plugin-react'

export default defineConfig({
  plugins: [react()],
  build: {
    outDir: 'dist',
    sourcemap: false,
    minify: 'terser',
    terserOptions: {
      compress: {
        drop_console: true,      // 删除所有 console.log
        drop_debugger: true,
        passes: 3,               // 三轮压缩，更彻底
        unsafe_math: true,
        pure_getters: true,
      },
      mangle: {
        toplevel: true,          // 顶层变量名也混淆
        properties: {
          // 不混淆属性名，否则会破坏 React/Recharts
          regex: /^_/            // 只混淆以_开头的私有属性
        }
      },
      format: {
        comments: false,         // 删除所有注释
        ascii_only: true,        // 非ASCII字符转义（中文→\uXXXX）
      }
    },
    rollupOptions: {
      output: {
        // 文件名加hash，防止缓存旧版本
        entryFileNames: 'assets/[name]-[hash].js',
        chunkFileNames: 'assets/[name]-[hash].js',
        assetFileNames: 'assets/[name]-[hash].[ext]',
      }
    }
  }
})
