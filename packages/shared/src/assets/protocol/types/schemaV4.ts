import type { DeckSlotId, ModuleModel } from '../../../types'
import type {
  ProtocolFile as V3ProtocolFile,
  AspDispAirgapParams,
  BlowoutParams,
  TouchTipParams,
  PipetteAccessParams,
  MoveToSlotParams,
  DelayParams,
} from './schemaV3'

export type { BlowoutParams, FilePipette, FileLabware } from './schemaV3'

export interface FileModule {
  slot: DeckSlotId
  model: ModuleModel
}

export interface EngageMagnetParams {
  module: string
  engageHeight: number
}

export interface TemperatureParams {
  module: string
  temperature: number
}

export interface AtomicProfileStep {
  holdTime: number
  temperature: number
}

export interface TCProfileParams {
  module: string
  profile: AtomicProfileStep[]
  volume: number
}

export interface ModuleOnlyParams {
  module: string
}

export interface ThermocyclerSetBlockTemperatureArgs {
  module: string
  temperature: number
  volume?: number
}

export type Command =
  | {
      command: 'aspirate' | 'dispense' | 'airGap'
      params: AspDispAirgapParams
    }
  | {
      command: 'blowout'
      params: BlowoutParams
    }
  | {
      command: 'touchTip'
      params: TouchTipParams
    }
  | {
      command: 'pickUpTip' | 'dropTip'
      params: PipetteAccessParams
    }
  | {
      command: 'moveToSlot'
      params: MoveToSlotParams
    }
  | {
      command: 'delay'
      params: DelayParams
    }
  | {
      command: 'magneticModule/engageMagnet'
      params: EngageMagnetParams
    }
  | {
      command: 'magneticModule/disengageMagnet'
      params: ModuleOnlyParams
    }
  | {
      command: 'temperatureModule/setTargetTemperature'
      params: TemperatureParams
    }
  | {
      command: 'temperatureModule/deactivate'
      params: ModuleOnlyParams
    }
  | {
      command: 'temperatureModule/awaitTemperature'
      params: TemperatureParams
    }
  | {
      command: 'thermocycler/setTargetBlockTemperature'
      params: ThermocyclerSetBlockTemperatureArgs
    }
  | {
      command: 'thermocycler/setTargetLidTemperature'
      params: TemperatureParams
    }
  | {
      command: 'thermocycler/awaitBlockTemperature'
      params: TemperatureParams
    }
  | {
      command: 'thermocycler/awaitLidTemperature'
      params: TemperatureParams
    }
  | {
      command: 'thermocycler/openLid'
      params: ModuleOnlyParams
    }
  | {
      command: 'thermocycler/closeLid'
      params: ModuleOnlyParams
    }
  | {
      command: 'thermocycler/deactivateBlock'
      params: ModuleOnlyParams
    }
  | {
      command: 'thermocycler/deactivateLid'
      params: ModuleOnlyParams
    }
  | {
      command: 'thermocycler/runProfile'
      params: TCProfileParams
    }
  | {
      command: 'thermocycler/awaitProfileComplete'
      params: ModuleOnlyParams
    }

// NOTE: must be kept in sync with '../schemas/4.json'
export interface ProtocolFile<DesignerApplicationData>
  extends Omit<
    V3ProtocolFile<DesignerApplicationData>,
    'schemaVersion' | 'commands'
  > {
  $otSharedSchema: '#/protocol/schemas/4'
  schemaVersion: 4
  modules: Record<string, FileModule>
  commands: Command[]
}
