import { createRouter, createWebHistory } from 'vue-router'
import HomeView from '@/views/HomeView.vue'
import Main from '@/pages/main/index.vue'
import Cookies from 'js-cookie'

const router = createRouter({
  history: createWebHistory(import.meta.env.BASE_URL),
  routes: [
    {
      path: '/',
      name: 'Home',
      component: Main,
    },
    {
      path: '/login',
      name: 'Login',
      component: () => import('@/pages/login/index.vue')
    },
    {
      path: '/about',
      name: 'about',
      // route level code-splitting
      // this generates a separate chunk (About.[hash].js) for this route
      // which is lazy-loaded when the route is visited.
      component: () => import('../views/AboutView.vue')
    }
  ]
})

router.beforeEach(async (to, from) => {
  const token = Cookies.get('x-token');
  // console.log('====', token)
  if (!token && to.name !== 'Login') {
    return { name: 'Login' }
  }
  return true;
})

export default router
