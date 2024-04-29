#Code block description: This script block is suitable for the make DNB steps in a two-step experiment. This script prepares an Opentrons OT-2 robot for a lab experiment by setting up a temperature module and a thermocycler module with associated labware. It also configures single and multi-channel pipettes with tip racks, ensuring the robot is ready for precise liquid handling tasks in a molecular biology protocol.

from opentrons.types import Point
metadata = {
    'protocolName': 'mtb',
    'author': 'MGI-X',
    'apiLevel': '2.14'
}

def run(ctx):
    # Module definition
    temp_mod = ctx.load_module('temperature module gen2', '6')
    temp_rack=temp_mod.load_labware('biorad_96_wellplate_200ul_pcr') 
    tc_mod = ctx.load_module(module_name='thermocyclerModuleV1')
    tc_rack = tc_mod.load_labware(name='biorad_96_wellplate_200ul_pcr')
    # pipette and tiprack
    tiprack = [ctx.load_labware('opentrons_96_filtertiprack_20ul', slot) for slot in ['3']]
    pip20_single = ctx.load_instrument('p20_single_gen2', 'right', tip_racks=[*tiprack])
    pip20_multi = ctx.load_instrument('p20_multi_gen2', 'left', tip_racks=[*tiprack])