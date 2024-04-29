#Code block description: This script block is suitable for the first round of purification steps in a two-step experiment. This script section handles the loading of reagents into specific wells of labware on an Opentrons OT-2 robot. TE buffer is loaded into well A1 of the temperature-controlled rack, PCR product into well A1 of the magnetic rack, ethanol into well A1 of a deep well plate, and beads into well B1 of the temperature-controlled rack. The volumes for each liquid are specified, ensuring precise reagent handling for subsequent experimental steps.
from opentrons.types import Point
metadata = {
    'protocolName': 'mtb',
    'author': 'MGI-X',
    'apiLevel': '2.14'
}

def run(ctx):
    #Reagent loading
    temp_rack["A1"].load_liquid(liquid=TE_buffer, volume = 13)
    mag_rack["A1"].load_liquid(liquid={round}st_pcr_product, volume = 25)
    plate_1["A1"].load_liquid(liquid=Ethonal, volume = 480)
    temp_rack["B1"].load_liquid(liquid=Beads, volume = 50)