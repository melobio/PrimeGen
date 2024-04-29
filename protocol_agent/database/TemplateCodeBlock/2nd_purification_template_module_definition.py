#Code block description: This script block is suitable for the second round of purification steps in a two-step experiment. This script sets up an Opentrons OT-2 robot for a laboratory protocol by configuring modules for temperature control and magnetic separation. It loads specific labware for reagents and samples, and prepares pipetting instruments with tip racks for precise liquid handling. Essential for automated biochemical assays.

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
    mag_mod = ctx.load_module('magnetic module gen2', '3')
    mag_rack = mag_mod.load_labware('biorad_96_wellplate_200ul_pcr')
    # Plates
    plate_1 = ctx.load_labware('nest_96_wellplate_2ml_deep', '1', 'reg1')
    # pipette and tiprack
    tiprack = [ctx.load_labware('opentrons_96_filtertiprack_20ul', slot) for slot in ['2']] 
    pip20_single = ctx.load_instrument('p20_single_gen2', 'right', tip_racks=[*tiprack]) 
    pip20_multi = ctx.load_instrument('p20_multi_gen2', 'left', tip_racks=[*tiprack]) 