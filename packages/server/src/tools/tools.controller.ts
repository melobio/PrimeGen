import { Controller, Get } from '@nestjs/common';
import { ToolsService } from './tools.service';
import { XResponse } from '../common/response';

@Controller('tools/template')
export class ToolsController {
  constructor(private readonly toolsService: ToolsService) {}

  @Get('agents')
  async getAgentTemplates() {
    return XResponse.Success(await this.toolsService.getAgentTemplates());
  }
  @Get('llms')
  async getLLMTemplates() {
    return XResponse.Success(await this.toolsService.getLLMTemplates());
  }
  @Get('devices')
  async getDeviceTemplates() {
    return XResponse.Success(await this.toolsService.getDeviceTemplates());
  }
  @Get('prompts')
  async getPromptTemplates() {
    return XResponse.Success(await this.toolsService.getPromptTemplates());
  }
  @Get('inputs')
  async getInputTemplates() {
    return XResponse.Success(await this.toolsService.getInputTemplates());
  }
  @Get('tools')
  async getToolTemplates() {
    return XResponse.Success(await this.toolsService.getToolTemplates());
  }
}
