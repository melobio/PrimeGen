#Code block description: This script block is suitable for the second round of PCR reaction steps in a two-step experiment. This script sets up an Opentrons OT-2 robot for a PCR experiment. It initializes a temperature module, a thermocycler, and loads PCR plates. It also calculates the number of tip racks needed and equips the robot with single and multi-channel pipettes, preparing it for precise liquid handling.

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
    # Plates
    plate_1 = ctx.load_labware('biorad_96_wellplate_200ul_pcr', '1', 'reg1') 
    # The number of tipracks needs to be calculated
    tiprack = [ctx.load_labware('opentrons_96_filtertiprack_20ul', slot) for slot in 2nd_Tip_rack_number] 
    pip20_single = ctx.load_instrument('p20_single_gen2', 'right', tip_racks=[*tiprack]) 
    pip20_multi = ctx.load_instrument('p20_multi_gen2', 'left', tip_racks=[*tiprack]) 
