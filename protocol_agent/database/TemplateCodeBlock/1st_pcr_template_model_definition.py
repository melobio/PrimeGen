#Code block description: This script block is suitable for the first round of PCR reaction steps in a two-step experiment. This script configures the hardware setup for an Opentrons OT-2 robot, preparing it for a molecular biology experiment. It loads a second-generation temperature module and a thermocycler module, each with a Bio-Rad 96-well plate for PCR reactions. It also loads another 96-well plate on the robot's deck. Tip racks are dynamically loaded based on specified slots, and two types of pipettes, a single-channel and a multi-channel, are equipped, each associated with the loaded tip racks. This setup is foundational for subsequent liquid handling and temperature control steps in the protocol.

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
    plate_1 = ctx.load_labware('biorad_96_wellplate_200ul_pcr', '1', 'reg1') 
    tiprack = [ctx.load_labware('opentrons_96_filtertiprack_20ul', slot) for slot in 1st_Tip_rack_number] 
    pip20_single = ctx.load_instrument('p20_single_gen2', 'right', tip_racks=[*tiprack]) 
    pip20_multi = ctx.load_instrument('p20_multi_gen2', 'left', tip_racks=[*tiprack]) 