import { DynamicModule, Module } from '@nestjs/common';
import { RedisClientOptions } from 'redis';
import { CacheService } from './cache.service';
import { StoreConfig } from 'cache-manager';

export type CacheOptions = RedisClientOptions & StoreConfig & {
  isGlobal: boolean;
  ttl?: number;
};

export const CACHE_MODULE = 'CACHE_MODULE';

@Module({})
export class CacheModule {
  static register(options: CacheOptions): DynamicModule {
    return {
      module: CacheModule,
      providers: [
        {
          provide: 'CONFIG_OPTIONS',
          useValue: options,
        },
        {
          provide: CACHE_MODULE,
          useClass: CacheService,
        },
      ],
      exports: [CACHE_MODULE],
      global: options.isGlobal,
    };
  }
}
