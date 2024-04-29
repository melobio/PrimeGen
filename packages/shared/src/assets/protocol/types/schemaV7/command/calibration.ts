import type { CommonCommandRunTimeInfo, CommonCommandCreateInfo } from './index'
import type {
  PipetteMount,
  GantryMount,
  LabwareOffset,
  Coordinates,
} from '../../../../../types'
// TODO (sb 10/26/22): Separate out calibration commands from protocol schema in RAUT-272
export interface CalibratePipetteCreateCommand extends CommonCommandCreateInfo {
  commandType: 'calibration/calibratePipette'
  params: CalibratePipetteParams
}
export interface CalibrateGripperCreateCommand extends CommonCommandCreateInfo {
  commandType: 'calibration/calibrateGripper'
  params: CalibrateGripperParams
}
export interface CalibrateModuleCreateCommand extends CommonCommandCreateInfo {
  commandType: 'calibration/calibrateModule'
  params: CalibrateModuleParams
}
export interface MoveToMaintenancePositionCreateCommand
  extends CommonCommandCreateInfo {
  commandType: 'calibration/moveToMaintenancePosition'
  params: MoveToMaintenancePositionParams
}

export interface CalibratePipetteRunTimeCommand
  extends CommonCommandRunTimeInfo,
    CalibratePipetteCreateCommand {
  result?: CalibratePipetteResult
}
export interface CalibrateGripperRunTimeCommand
  extends CommonCommandRunTimeInfo,
    CalibrateGripperCreateCommand {
  result?: CalibrateGripperResult
}
export interface CalibrateModuleRunTimeCommand
  extends CommonCommandRunTimeInfo,
    CalibrateModuleCreateCommand {
  result?: CalibrateModuleResult
}
export interface MoveToMaintenancePositionRunTimeCommand
  extends CommonCommandRunTimeInfo,
    MoveToMaintenancePositionCreateCommand {
  result?: {}
}

export type CalibrationRunTimeCommand =
  | CalibratePipetteRunTimeCommand
  | CalibrateGripperRunTimeCommand
  | CalibrateModuleRunTimeCommand
  | MoveToMaintenancePositionRunTimeCommand

export type CalibrationCreateCommand =
  | CalibratePipetteCreateCommand
  | CalibrateGripperCreateCommand
  | CalibrateModuleCreateCommand
  | MoveToMaintenancePositionCreateCommand

interface CalibratePipetteParams {
  mount: PipetteMount
}
interface CalibrateGripperParams {
  jaw: 'front' | 'rear'
  otherJawOffset?: Coordinates
}
interface CalibrateModuleParams {
  moduleId: string
  labwareId: string
  mount: PipetteMount
}
interface CalibratePipetteResult {
  pipetteOffset: LabwareOffset
}
interface CalibrateGripperResult {
  jawOffset: Coordinates
}
interface MoveToMaintenancePositionParams {
  mount: GantryMount
  maintenancePosition?: 'attachPlate' | 'attachInstrument'
}
interface CalibrateModuleResult {
  moduleOffset: Coordinates
}
