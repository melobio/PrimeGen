<template>
  <div class="ot2-commands">
    <div class="commands-item"
         :key="command.id"
         v-for="(command, index) in commandRun">
      <!--{{ command.commandType }}-->
      <div>
        <span
          class="command-index">{{ index + 1 }}.</span>
        <CommandText :command="command"/>
      </div>
      <JetFaultCheckResult
        v-if="commandIndexToFailCheckResult[index]"
        v-bind="commandIndexToFailCheckResult[index]"
        :img-width="300"
      />
    </div>
  </div>
</template>

<script setup lang='ts'>
import type {PropType} from "vue";
import JetFaultCheckResult from "@/pages/main/chat/chat-detail/JetFaultCheckResult.vue";
import * as _ from 'lodash'
import {computed, h} from 'vue'
import { useTranslation } from "i18next-vue";
const { t } = useTranslation('protocol_command_text');
import type {
  AspirateRunTimeCommand,
  BlowoutRunTimeCommand,
  CompletedProtocolAnalysis,
  DispenseRunTimeCommand,
  DropTipRunTimeCommand, HeaterShakerSetTargetTemperatureCreateCommand,
  MoveToWellRunTimeCommand,
  PickUpTipRunTimeCommand,
  RunTimeCommand, TCSetTargetBlockTemperatureCreateCommand, TCSetTargetLidTemperatureCreateCommand,
  TemperatureModuleAwaitTemperatureCreateCommand,
  TemperatureModuleSetTargetTemperatureCreateCommand
} from "@xpcr/shared";
import {getLoadedLabware} from "@/pages/main/chat/chat-detail/accessors";
import {
  getLabwareDefinitionsFromCommands,
  getLabwareDisplayLocation
} from "@/pages/main/chat/chat-detail/getLabwareDisplayLocation";
import {
  getModuleModel
} from "@/pages/main/chat/chat-detail/getModuleModel";
import {
  getPipetteNameOnMount
} from "@/pages/main/chat/chat-detail/getPipetteNameOnMount";
import {
  getModuleDisplayName,
  getLabwareDefURI,
  getModuleType,
  getOccludedSlotCountForModule,
  getPipetteNameSpecs,
  OT2_STANDARD_MODEL,
} from "@xpcr/shared/src/index";
import {
  getModuleDisplayLocation
} from "@/pages/main/chat/chat-detail/getModuleDisplayLocation";
import {
  getFinalLabwareLocation
} from "@/pages/main/chat/chat-detail/getFinalLabwareLocation";
import type {
  LabwareLocation,
  LoadLabwareRunTimeCommand,
  MoveLabwareRunTimeCommand,
} from "@xpcr/shared";
import { getLabwareName } from "@/pages/main/chat/chat-detail/getLabwareName";
import { getLiquidDisplayName } from '@/pages/main/chat/chat-detail/getLiquidDisplayName'

const props = defineProps({
  protocolAnalysis: {
    type: Object as PropType<CompletedProtocolAnalysis>,
    // type: Object,
    required: true
  },
  currentCommandIndex: {
    type: Number,
    default: -1,
  },
  commandFailCheckResult: {
    type: Array,
    default: () => [],
  },
});

// const commandIndexToFailCheckResult: { [commandIndex: number]: any } = {}
// props.commandFailCheckResult.forEach((a: any) => {
//   commandIndexToFailCheckResult[a.commandIndex] = a
// })
const commandIndexToFailCheckResult = computed(() => {
  const result: { [commandIndex: number]: any } = {}
  props.commandFailCheckResult.forEach((a: any) => {
    result[a.commandIndex] = a
  })
  return result
})

const commandRun = computed(() => {
  return _.slice(props.protocolAnalysis.commands, 0, props.currentCommandIndex + 1);
});

const SIMPLE_TRANSLATION_KEY_BY_COMMAND_TYPE: {
  [commandType in RunTimeCommand['commandType']]?: string
} = {
  home: 'home_gantry',
  savePosition: 'save_position',
  touchTip: 'touch_tip',
  'magneticModule/engage': 'engaging_magnetic_module',
  'magneticModule/disengage': 'disengaging_magnetic_module',
  'temperatureModule/deactivate': 'deactivate_temperature_module',
  'thermocycler/waitForBlockTemperature': 'waiting_for_tc_block_to_reach',
  'thermocycler/waitForLidTemperature': 'waiting_for_tc_lid_to_reach',
  'thermocycler/openLid': 'opening_tc_lid',
  'thermocycler/closeLid': 'closing_tc_lid',
  'thermocycler/deactivateBlock': 'deactivating_tc_block',
  'thermocycler/deactivateLid': 'deactivating_tc_lid',
  'thermocycler/awaitProfileComplete': 'tc_awaiting_for_duration',
  'heaterShaker/deactivateHeater': 'deactivating_hs_heater',
  'heaterShaker/openLabwareLatch': 'unlatching_hs_latch',
  'heaterShaker/closeLabwareLatch': 'latching_hs_latch',
  'heaterShaker/deactivateShaker': 'deactivate_hs_shake',
  'heaterShaker/waitForTemperature': 'waiting_for_hs_to_reach',
}

const robotSideAnalysis = props.protocolAnalysis;

const CommandText = (props: { command: RunTimeCommand }) => {
  const { command } = props;
  const { t } = useTranslation('protocol_command_text')
  switch (command.commandType) {
    case 'aspirate':
    case 'dispense':
    case 'blowout':
    case 'moveToWell':
    case 'dropTip':
    case 'pickUpTip': {
      return h('span', {}, PipettingCommandText({ command })); // PipettingCommandText({ command })
    }
    case 'loadLabware':
    case 'loadPipette':
    case 'loadModule':
    case 'loadLiquid': {
      return h('span', {}, LoadCommandText({ command}));
    }
    case 'temperatureModule/setTargetTemperature':
    case 'temperatureModule/waitForTemperature':
    case 'thermocycler/setTargetBlockTemperature':
    case 'thermocycler/setTargetLidTemperature':
    case 'heaterShaker/setTargetTemperature': {
      return h('span', {}, TemperatureCommandText({ command }));
    }
    case 'thermocycler/runProfile': {
      const { profile } = command.params
      const steps = profile.map(
        ({ holdSeconds, celsius }: { holdSeconds: number; celsius: number }) =>
          t('tc_run_profile_steps', { celsius: celsius, seconds: holdSeconds })
      )
      const startingProfile = t('tc_starting_profile', {
        repetitions: Object.keys(steps).length,
      });
      const children = []
      children.push(h('div', {}, startingProfile))
      steps.forEach((step, index) => {
        children.push(h('div', {}, step))
      })
      return h('div', {}, children);
    }
    case 'heaterShaker/setAndWaitForShakeSpeed': {
      const { rpm } = command.params
      return h('span', {}, t('set_and_await_hs_shake', { rpm }));
    }
    case 'moveToSlot': {
      const { slotName } = command.params
      return h('span', {}, t('move_to_slot', { slot_name: slotName }));
    }
    case 'moveRelative': {
      const { axis, distance } = command.params
      return h('span', {}, t('move_relative', { axis, distance }));
    }
    case 'moveToCoordinates': {
      const { coordinates } = command.params
      // t('move_to_coordinates', coordinates)
      return h('span', {}, t('move_to_coordinates', { coordinates }));
    }
    case 'moveLabware': {
      return h('span', {}, MoveLabwareCommandText({ command }));
    }
    case 'configureForVolume': {
      const { volume, pipetteId } = command.params
      const pipetteName = robotSideAnalysis.pipettes.find(
        pip => pip.id === pipetteId
      )?.pipetteName

      return h('span', {}, t('configure_for_volume', {
        volume,
        pipette:
          pipetteName != null
            ? getPipetteNameSpecs(pipetteName)?.displayName
            : '',
      }));
    }
    case 'touchTip':
    case 'home':
    case 'savePosition':
    case 'magneticModule/engage':
    case 'magneticModule/disengage':
    case 'temperatureModule/deactivate':
    case 'thermocycler/waitForBlockTemperature':
    case 'thermocycler/waitForLidTemperature':
    case 'thermocycler/openLid':
    case 'thermocycler/closeLid':
    case 'thermocycler/deactivateBlock':
    case 'thermocycler/deactivateLid':
    case 'thermocycler/awaitProfileComplete':
    case 'heaterShaker/deactivateHeater':
    case 'heaterShaker/openLabwareLatch':
    case 'heaterShaker/closeLabwareLatch':
    case 'heaterShaker/deactivateShaker':
    case 'heaterShaker/waitForTemperature': {
      const simpleTKey =
        SIMPLE_TRANSLATION_KEY_BY_COMMAND_TYPE[command.commandType]
      return h('span', {}, simpleTKey != null ? t(simpleTKey) : '');
    }
    case 'waitForDuration': {
      const { seconds, message } = command.params
      // t('wait_for_duration', { seconds, message })
      return h('span', {}, t('waiting_for_duration', { seconds, message }));
    }
    case 'pause': // legacy pause command
    case 'waitForResume': {
      return h('span', {}, command.params?.message && command.params.message !== ''
        ? command.params.message
        : t('wait_for_resume'));
    }
    case 'delay': {
      const { message = '' } = command.params
      if ('waitForResume' in command.params) {
        return h('span', {}, command.params?.message && command.params.message !== ''
          ? command.params.message
          : t('wait_for_resume'));
      } else {
        return h('span', {}, t('wait_for_duration', {
          seconds: command.params.seconds,
          message,
        }));
      }
    }
    case 'custom': {
      const { legacyCommandText } = command.params ?? {}
      const sanitizedCommandText =
        typeof legacyCommandText === 'object'
          ? JSON.stringify(legacyCommandText)
          : String(legacyCommandText)
      return h('span', {}, legacyCommandText != null
        ? sanitizedCommandText
        : `${command.commandType}: ${JSON.stringify(command.params)}`);
    }
    default: {
      console.warn(
        'CommandText encountered a command with an unrecognized commandType: ',
        command
      )
      return h('span', {}, JSON.stringify(command));
    }
  }
}
CommandText.props = ['command'];

type PipettingRunTimeCommmand =
  | AspirateRunTimeCommand
  | DispenseRunTimeCommand
  | BlowoutRunTimeCommand
  | MoveToWellRunTimeCommand
  | DropTipRunTimeCommand
  | PickUpTipRunTimeCommand

const PipettingCommandText = (props: { command: PipettingRunTimeCommmand }): string => {
  const { command } = props;
  const { labwareId, wellName } = command.params
  const labwareLocation = getLoadedLabware(robotSideAnalysis, labwareId)
    ?.location
  const displayLocation =
    labwareLocation
      ? getLabwareDisplayLocation(robotSideAnalysis, labwareLocation, t)
      : ''


  switch (command.commandType) {
    case 'aspirate': {
      const { volume, flowRate } = command.params
      return t('aspirate', {
        well_name: wellName,
        labware: getLabwareName(robotSideAnalysis, labwareId),
        labware_location: displayLocation,
        volume: volume,
        flow_rate: flowRate,
      })
    }
    case 'dispense': {
      const { volume, flowRate, pushOut } = command.params
      return pushOut
        ? t('dispense_push_out', {
          well_name: wellName,
          labware: getLabwareName(robotSideAnalysis, labwareId),
          labware_location: displayLocation,
          volume: volume,
          flow_rate: flowRate,
          push_out_volume: pushOut,
        })
        : t('dispense', {
          well_name: wellName,
          labware: getLabwareName(robotSideAnalysis, labwareId),
          labware_location: displayLocation,
          volume: volume,
          flow_rate: flowRate,
        })
    }
    case 'blowout': {
      const { flowRate } = command.params
      return t('blowout', {
        well_name: wellName,
        labware: getLabwareName(robotSideAnalysis, labwareId),
        labware_location: displayLocation,
        flow_rate: flowRate,
      })
    }
    case 'moveToWell': {
      return t('move_to_well', {
        well_name: wellName,
        labware: getLabwareName(robotSideAnalysis, labwareId),
        labware_location: displayLocation,
      })
    }
    case 'dropTip': {
      const loadedLabware = getLoadedLabware(robotSideAnalysis, labwareId)
      const labwareDefinitions = getLabwareDefinitionsFromCommands(
        robotSideAnalysis.commands
      )
      const labwareDef = labwareDefinitions.find(
        lw => getLabwareDefURI(lw) === loadedLabware?.definitionUri
      )
      return labwareDef?.parameters.isTiprack
        ? t('return_tip', {
          well_name: wellName,
          labware: getLabwareName(robotSideAnalysis, labwareId),
          labware_location: displayLocation,
        })
        : t('drop_tip', {
          well_name: wellName,
          labware: getLabwareName(robotSideAnalysis, labwareId),
        })
    }
    case 'pickUpTip': {
      return t('pickup_tip', {
        well_name: wellName,
        labware: getLabwareName(robotSideAnalysis, labwareId),
        labware_location: displayLocation,
      })
    }
    default: {
      console.warn(
        'PipettingCommandText encountered a command with an unrecognized commandType: ',
        command
      )
      return ''
    }
  }
}

PipettingCommandText.props = ['command']
const LoadCommandText = (props: { command: RunTimeCommand; }): string => {
  const { command } = props
  const { t } = useTranslation('run_details')
  switch (command.commandType) {
    case 'loadPipette': {
      const pipetteModel = getPipetteNameOnMount(
        robotSideAnalysis,
        command.params.mount
      )
      return t('load_pipette_protocol_setup', {
        pipette_name:
          pipetteModel != null
            ? getPipetteNameSpecs(pipetteModel)?.displayName ?? ''
            : '',
        mount_name: command.params.mount === 'left' ? t('left') : t('right'),
      })
    }
    case 'loadModule': {
      const occludedSlotCount = getOccludedSlotCountForModule(
        getModuleType(command.params.model),
        robotSideAnalysis.robotType ?? OT2_STANDARD_MODEL
      )
      return t('load_module_protocol_setup', {
        count: occludedSlotCount,
        module: getModuleDisplayName(command.params.model),
        slot_name: command.params.location.slotName,
      })
    }
    case 'loadLabware': {
      if (
        command.params.location !== 'offDeck' &&
        'moduleId' in command.params.location
      ) {
        const moduleModel = getModuleModel(
          robotSideAnalysis,
          command.params.location.moduleId
        )
        const moduleName =
          moduleModel != null ? getModuleDisplayName(moduleModel) : ''

        return t('load_labware_info_protocol_setup', {
          count:
            moduleModel != null
              ? getOccludedSlotCountForModule(
                getModuleType(moduleModel),
                robotSideAnalysis.robotType ?? OT2_STANDARD_MODEL
              )
              : 1,
          labware: command.result?.definition.metadata.displayName,
          slot_name: getModuleDisplayLocation(
            robotSideAnalysis,
            command.params.location.moduleId
          ),
          module_name: moduleName,
        })
      } else if (
        command.params.location !== 'offDeck' &&
        'labwareId' in command.params.location
      ) {
        const labwareId = command.params.location.labwareId
        const labwareName = command.result?.definition.metadata.displayName
        const matchingAdapter = robotSideAnalysis.commands.find(
          (command): command is LoadLabwareRunTimeCommand =>
            command.commandType === 'loadLabware' &&
            command.result?.labwareId === labwareId
        )
        const adapterName =
          matchingAdapter?.result?.definition.metadata.displayName
        const adapterLoc = matchingAdapter?.params.location
        if (adapterLoc === 'offDeck') {
          return t('load_labware_info_protocol_setup_adapter_off_deck', {
            labware: labwareName,
            adapter_name: adapterName,
          })
        } else if (adapterLoc != null && 'slotName' in adapterLoc) {
          return t('load_labware_info_protocol_setup_adapter', {
            labware: labwareName,
            adapter_name: adapterName,
            slot_name: adapterLoc?.slotName,
          })
        } else if (adapterLoc != null && 'moduleId' in adapterLoc) {
          const moduleModel = getModuleModel(
            robotSideAnalysis,
            adapterLoc?.moduleId ?? ''
          )
          const moduleName =
            moduleModel != null ? getModuleDisplayName(moduleModel) : ''
          return t('load_labware_info_protocol_setup_adapter_module', {
            labware: labwareName,
            adapter_name: adapterName,
            module_name: moduleName,
            slot_name: getModuleDisplayLocation(
              robotSideAnalysis,
              adapterLoc?.moduleId ?? ''
            ),
          })
        } else {
          //  shouldn't reach here, adapter shouldn't have location  type labwareId
          return ''
        }
      } else {
        const labware = command.result?.definition.metadata.displayName
        return command.params.location === 'offDeck'
          ? t('load_labware_info_protocol_setup_off_deck', { labware })
          : t('load_labware_info_protocol_setup_no_module', {
            labware,
            slot_name: command.params.location?.slotName,
          })
      }
    }
    case 'loadLiquid': {
      const { liquidId, labwareId } = command.params
      return t('load_liquids_info_protocol_setup', {
        liquid: getLiquidDisplayName(robotSideAnalysis, liquidId),
        labware: getLabwareName(robotSideAnalysis, labwareId),
      })
    }
    default: {
      console.warn(
        'LoadCommandText encountered a command with an unrecognized commandType: ',
        command
      )
      return ''
    }
  }
}
LoadCommandText.props = ['command']

type TemperatureCreateCommand =
  | TemperatureModuleSetTargetTemperatureCreateCommand
  | TemperatureModuleAwaitTemperatureCreateCommand
  | TCSetTargetBlockTemperatureCreateCommand
  | TCSetTargetLidTemperatureCreateCommand
  | HeaterShakerSetTargetTemperatureCreateCommand

const T_KEYS_BY_COMMAND_TYPE: {
  [commandType in TemperatureCreateCommand['commandType']]: string
} = {
  'temperatureModule/setTargetTemperature': 'setting_temperature_module_temp',
  'temperatureModule/waitForTemperature': 'waiting_to_reach_temp_module',
  'thermocycler/setTargetBlockTemperature': 'setting_thermocycler_block_temp',
  'thermocycler/setTargetLidTemperature': 'setting_thermocycler_lid_temp',
  'heaterShaker/setTargetTemperature': 'setting_hs_temp',
}

const TemperatureCommandText = (props: {command: TemperatureCreateCommand}) => {
  const { command } = props
  const { t } = useTranslation('protocol_command_text')

  return t(T_KEYS_BY_COMMAND_TYPE[command.commandType] as any, {
    temp:
      command.params?.celsius != null
        ? t('degrees_c', { temp: command.params.celsius })
        : t('target_temperature'),
    hold_time_seconds:
      'holdTimeSeconds' in command.params
        ? command.params.holdTimeSeconds ?? '0'
        : '0',
  })
}
TemperatureCommandText.props = ['command']

const MoveLabwareCommandText = (props: { command: MoveLabwareRunTimeCommand }): string => {
  const { t } = useTranslation('protocol_command_text')
  const { command } = props
  const { labwareId, newLocation, strategy } = command.params

  const allPreviousCommands = robotSideAnalysis.commands.slice(
    0,
    robotSideAnalysis.commands.findIndex(c => c.id === command.id)
  )
  const oldLocation = getFinalLabwareLocation(labwareId, allPreviousCommands)
  const newDisplayLocation = getLabwareDisplayLocation(
    robotSideAnalysis,
    newLocation,
    t
  )

  return strategy === 'usingGripper'
    ? t('move_labware_using_gripper', {
      labware: getLabwareName(robotSideAnalysis, labwareId),
      old_location:
        oldLocation != null
          ? getLabwareDisplayLocation(robotSideAnalysis, oldLocation, t)
          : '',
      new_location: newDisplayLocation,
    })
    : t('move_labware_manually', {
      labware: getLabwareName(robotSideAnalysis, labwareId),
      old_location:
        oldLocation != null
          ? getLabwareDisplayLocation(robotSideAnalysis, oldLocation, t)
          : '',
      new_location: newDisplayLocation,
    })
}

MoveLabwareCommandText.props = ['command'];
</script>

<style scoped lang='scss'>
.ot2-commands {
  display: flex;
  flex-direction: column;
  .commands-item {
    display: flex;
    flex-direction: column;
    align-items: start;
    margin-bottom: 5px;
    .command-index {
      color: white;
      font-weight: bold;
      font-size: 14px;
      margin-right: 5px;
    }
  }
}
</style>