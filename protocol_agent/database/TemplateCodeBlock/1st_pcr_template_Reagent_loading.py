#Code block description: This script block is suitable for the first round of PCR reaction steps in a two-step experiment. This script is for an Opentrons OT-2 lab automation robot, designed for automated reagent and sample loading. It involves loading specific PCR barcode reactants into a 96-well plate, PCR enzyme mixes and primer pools into a temporary rack, and DNA templates into a temperature-controlled rack, with each step specifying the type of liquid and volume.

from opentrons.types import Point
metadata = {
    'protocolName': 'mtb',
    'author': 'MGI-X',
    'apiLevel': '2.14'
}

def run(ctx):
    #Reagent loading
    for well in plate_1.wells()[:96]:
        well.load_liquid(liquid=Barcode_{PCR_round}, volume=30) 

    for well_name in {1st_PCR_Mix_well}：
         temp_rack[well_name].load_liquid(liquid={PCR_round}st_PCR_Enzyme_Mix, volume=1st_PCR_MIX_Volume_fault_tolerance_volume) 
    
    temp_rack[temp_rack.wells()[len(1st_PCR_Mix_well)+1].load_liquid(liquid=primer_pools, volume=Primer_Pool_well_Volume) 
    
    for well_name in DNA_TEMPLATE_WELL: 
        tc_rack.wells_by_name()[well_name].load_liquid(liquid=DNA_TEMPLATE, volume=2.5) 