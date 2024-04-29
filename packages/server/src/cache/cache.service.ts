import { Inject, Injectable } from '@nestjs/common';
import { RedisStore, redisStore } from 'cache-manager-redis-store';
import { CacheOptions } from './cache.module';

@Injectable()
export class CacheService {
  private store: RedisStore;
  constructor(@Inject('CONFIG_OPTIONS') private options: CacheOptions) {
    redisStore(options).then((store) => {
      this.store = store;
    });
  }

  async setToken(key: string, value: any, ttl?: number) {
    await this.store.set(
      `token:${key}`,
      value,
      { ttl: ttl ?? this.options.ttl },
      undefined,
    );
  }

  async getToken<T>(key: string): Promise<T | undefined> {
    // get: <T>(key: string) => Promise<T | undefined>
    return await this.store.get(`token:${key}`, undefined, undefined);
  }

  async setObject<T>(key : string, value : T, ttl?: number) {
    await this.store.set(key, value, { ttl: ttl ?? this.options.ttl },
      undefined,);
  }

  async getObject<T>(key: string): Promise<T | undefined> {
    return await this.store.get(key, undefined, undefined);
  }
}
