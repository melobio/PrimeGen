import { Logger } from '@nestjs/common';
import { DeviceInterface } from './device-interface';

export abstract class BaseDevices implements DeviceInterface {
  get logger(): Logger {
    return this.getLogger();
  }
  abstract getLogger(): Logger;

  async init() {
    this.logger.log(`init ${this.constructor.name} success`);
  }
}
