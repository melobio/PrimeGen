import { defineStore } from 'pinia'
import { ref } from 'vue'
import {
  experiments as getExperiments,
  addConWithExp,
} from '@/api';
import type { ExperimentInterface } from "@xpcr/common/src/experiment";

export const useExperimentsStore = defineStore('experiments', () => {
  const experiments = ref<ExperimentInterface[]>([]);
  const currentExperiment = ref<ExperimentInterface>();
  const activeNodeId = ref<string>('');

  async function getAllExperiments() {
    const { data } =  await getExperiments<ExperimentInterface[]>()
      experiments.value = data || [];
      // 默认显示第一个实验
      if (!currentExperiment.value && experiments.value.length > 0) {
        currentExperiment.value = experiments.value[0];
      }
  }

  async function handleAddConWithExp() {
    return await addConWithExp()
  }
  function setCurrentExperiment(experiment: ExperimentInterface) {
    currentExperiment.value = experiment;
  }
  function setActiveNodeId(nodeId: string) {
    activeNodeId.value = nodeId;
  }

  return {
    experiments,
    currentExperiment,
    setCurrentExperiment,
    activeNodeId,
    setActiveNodeId,
    getAllExperiments,
    handleAddConWithExp
  }
});