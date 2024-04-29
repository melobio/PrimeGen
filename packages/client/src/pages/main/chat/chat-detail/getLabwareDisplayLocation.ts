import {
  getLabwareDisplayName,
  getLabwareDefURI,
  getModuleDisplayName,
  getModuleType,
  getOccludedSlotCountForModule,
  OT2_STANDARD_MODEL,
} from '@xpcr/shared/src/index'
import type { LabwareLocation, LabwareDefinition2 } from '@xpcr/shared'
import { getModuleDisplayLocation } from './getModuleDisplayLocation'
import { getModuleModel } from './getModuleModel'
// import { getLabwareDefinitionsFromCommands } from '../../LabwarePositionCheck/utils/labware'
import type { CompletedProtocolAnalysis } from '@xpcr/shared'
import type { TFunction } from 'i18next'
import type {
  ProtocolAnalysisOutput,
  RunTimeCommand,
} from '@xpcr/shared';


export function getLabwareDefinitionsFromCommands(
  commands: RunTimeCommand[]
): LabwareDefinition2[] {
  return commands.reduce<LabwareDefinition2[]>((acc, command) => {
    const isLoadingNewDef =
      command.commandType === 'loadLabware' &&
      !acc.some(
        def =>
          command.result?.definition != null &&
          getLabwareDefURI(def) === getLabwareDefURI(command.result?.definition)
      )

    return isLoadingNewDef && command.result?.definition != null
      ? [...acc, command.result?.definition]
      : acc
  }, [])
}

export function getLabwareDisplayLocation(
  robotSideAnalysis: CompletedProtocolAnalysis,
  location: LabwareLocation,
  t: TFunction<'protocol_command_text'>,
  isOnDevice?: boolean
): string {
  if (location === 'offDeck') {
    return t('off_deck')
  } else if ('slotName' in location) {
    return isOnDevice
      ? location.slotName
      : t('slot', { slot_name: location.slotName })
  } else if ('moduleId' in location) {
    const moduleModel = getModuleModel(robotSideAnalysis, location.moduleId)
    if (moduleModel == null) {
      console.warn('labware is located on an unknown module model')
      return ''
    } else {
      const slotName = getModuleDisplayLocation(
        robotSideAnalysis,
        location.moduleId
      )
      return isOnDevice
        ? `${getModuleDisplayName(moduleModel)}, ${slotName}`
        : t('module_in_slot', {
            count: getOccludedSlotCountForModule(
              getModuleType(moduleModel),
              robotSideAnalysis.robotType ?? OT2_STANDARD_MODEL
            ),
            module: getModuleDisplayName(moduleModel),
            slot_name: slotName,
          })
    }
  } else if ('labwareId' in location) {
    const adapter = robotSideAnalysis.labware.find(
      lw => lw.id === location.labwareId
    )
    const allDefs = getLabwareDefinitionsFromCommands(
      robotSideAnalysis.commands
    )
    const adapterDef = allDefs.find(
      def => getLabwareDefURI(def) === adapter?.definitionUri
    )
    const adapterDisplayName =
      adapterDef != null ? getLabwareDisplayName(adapterDef) : ''

    if (adapter == null) {
      console.warn('labware is located on an unknown adapter')
      return ''
    } else if (adapter.location === 'offDeck') {
      return t('off_deck')
    } else if ('slotName' in adapter.location) {
      return t('adapter_in_slot', {
        adapter: adapterDisplayName,
        slot: adapter.location.slotName,
      })
    } else if ('moduleId' in adapter.location) {
      const moduleIdUnderAdapter = adapter.location.moduleId
      const moduleModel = robotSideAnalysis.modules.find(
        module => module.id === moduleIdUnderAdapter
      )?.model
      if (moduleModel == null) {
        console.warn('labware is located on an adapter on an unknown module')
        return ''
      }
      const slotName = getModuleDisplayLocation(
        robotSideAnalysis,
        adapter.location.moduleId
      )
      return t('adapter_in_mod_in_slot', {
        count: getOccludedSlotCountForModule(
          getModuleType(moduleModel),
          robotSideAnalysis.robotType ?? OT2_STANDARD_MODEL
        ),
        module: getModuleDisplayName(moduleModel),
        adapter: adapterDisplayName,
        slot: slotName,
      })
    } else {
      console.warn(
        'display location on adapter could not be established: ',
        location
      )
      return ''
    }
  } else {
    console.warn('display location could not be established: ', location)
    return ''
  }
}
