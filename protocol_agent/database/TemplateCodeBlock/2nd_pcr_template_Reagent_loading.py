#Code block description: This script block is suitable for the second round of PCR reaction steps in a two-step experiment. This script involves loading reagents into specified wells of a plate setup on an Opentrons OT-2 robot. It automates the dispensing of barcode liquids, PCR enzyme mixes, UR solution, and DNA template into respective wells, using predefined volumes to ensure precise reaction setups for PCR analysis.

from opentrons.types import Point
metadata = {
    'protocolName': 'mtb',
    'author': 'MGI-X',
    'apiLevel': '2.14'
}

def run(ctx):
    #Reagent loading
    for well in plate_1.wells()[:48]:
        well.load_liquid(liquid=Barcode_{PCR_round}, volume=30) 
    for well_name in ["A1"]：
         temp_rack[well_name].load_liquid(liquid={PCR_round}st_PCR_Enzyme_Mix, volume=2nd_PCR_MIX_Volume_fault_tolerance_volume) 
    temp_rack[temp_rack.wells()[2].load_liquid(liquid=UR, volume=U-R_well_Volume) 
    for well_name in DNA_TEMPLATE_WELL: 
        tc_rack.wells_by_name()[well_name].load_liquid(liquid=DNA_TEMPLATE, volume=2.5) 
