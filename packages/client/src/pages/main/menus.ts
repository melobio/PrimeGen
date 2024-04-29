import IconChat from '@/assets/menu/chat.png'
import IconExperiment from '@/assets/menu/experiment.png'
import IconRobot from '@/assets/menu/robot.png'
import IconLabware from '@/assets/menu/labware.png'
import IconProtocol from '@/assets/menu/protocol.png'
import IconModule from '@/assets/menu/module.png'
import IconTools from '@/assets/menu/tools.png'
import IconKnowledge from '@/assets/menu/knowledge.png'
import IconModel from '@/assets/menu/model.png'
import IconHistory from '@/assets/menu/history.png'
export class Menus {
  constructor(
    readonly name: string,
    readonly icon: string,
    readonly group: string,
  ) {
  }
}

export const MENUS: Menus[] = [
  new Menus('Chats', IconChat, '1'),
  new Menus('Experiment', IconExperiment, '1'),
  new Menus('Robots', IconRobot, '2'),
  new Menus('Labware', IconLabware, '2'),
  new Menus('Protocol', IconProtocol, '2'),
  new Menus('Modules', IconModule, '2'),
  new Menus('Tools', IconTools, '3'),
  new Menus('Knowledge', IconKnowledge, '3'),
  new Menus('Model', IconModel, '3'),
  new Menus('History', IconHistory, '4'),
]