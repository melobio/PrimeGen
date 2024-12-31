import { Logger } from '@nestjs/common';
import { ConversationEntity } from '../../entities/conversation.entity';
import { DataSource } from 'typeorm';

export class PlanStep {
  constructor(
    public stepName: string,
    public purpose: string,
    public prompt: string,
    public state: string,
    public useAgents: string[],
    public checkWithContext: boolean,
    public purposeExample: string,
  ) {}
}
export class Planner {
  logger = new Logger(Planner.name);
  title = '';
  description = '';
  availableStates: string[] = [];
  steps: PlanStep[] = [];
  constructor(
    private cvs: ConversationEntity,
    private dataSource: DataSource,
  ) {}

  static parse(
    content: string,
    cvs: ConversationEntity,
    dataSource: DataSource,
  ): Planner {
    const obj = JSON.parse(content);
    const planner = new Planner(cvs, dataSource);
    planner.steps = obj.steps.map((step) => {
      return new PlanStep(
        step.stepName,
        step.purpose,
        step.prompt,
        step.state,
        step.useAgents,
        step.checkWithContext,
        step.example,
      );
    });
    planner.title = obj.title;
    planner.description = obj.description;
    planner.availableStates = obj.availableStates;
    if (planner.steps.length > 0) {
      const currentStepIndex = planner.steps.findIndex(
        (step) => step.stepName === cvs.currentStep,
      );
      planner.logger.debug('cvs.currentStep===>', cvs.currentStep);

      planner.logger.debug('currentStepIndex===>', currentStepIndex);
      if (currentStepIndex > 0) {
        planner.steps[currentStepIndex].state = 'running';
      } else {
        planner.steps[0].state = 'running';
        cvs.currentStep = planner.steps[0].stepName;
        dataSource.manager.update(ConversationEntity, cvs.id, cvs).then();
      }
    }
    return planner;
  }

  getCurrentStep(): PlanStep {
    return this.steps.find((step) => step.state === 'running');
  }

  getNextStep(): PlanStep {
    const currentStep = this.getCurrentStep();
    const currentIndex = this.steps.indexOf(currentStep);
    const nextStep = this.steps[currentIndex + 1];
    return nextStep;
  }

  goNextStep() {
    const currentStep = this.getCurrentStep();
    if (currentStep) {
      currentStep.state = 'finish';
      const currentIndex = this.steps.indexOf(currentStep);
      const nextStep = this.steps[currentIndex + 1];
      if (nextStep) {
        nextStep.state = 'running';
        this.cvs.currentStep = nextStep.stepName;
        this.dataSource.manager
          .update(ConversationEntity, this.cvs.id, this.cvs)
          .then();
      }
      return nextStep;
    }
  }

  setCurrentStep(currentStepName: string) {
    let currentStepIndex = -1;
    this.steps.forEach((step, index) => {
      if (step.stepName === currentStepName) {
        currentStepIndex = index;
      }
    });

    if (currentStepIndex > -1) {
      this.steps.forEach((step, index) => {
        if (index < currentStepIndex) {
          step.state = 'finish';
        } else if (index === currentStepIndex) {
          step.state = 'running';
        } else {
          step.state = 'idle';
        }
      });
      this.logger.log(`Current Step set to ${currentStepName}`);
      this.cvs.currentStep = currentStepName;
      this.dataSource.manager
        .update(ConversationEntity, this.cvs.id, this.cvs)
        .then();
    } else {
      this.logger.warn(`CurrentStepName ${currentStepName} not found `);
    }
  }
}
