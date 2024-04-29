#Code block description: This script block is suitable for the second round of purification steps in a two-step experiment. This script handles the loading of specific reagents into designated wells on an Opentrons OT-2 robot. TE buffer, PCR products, ethanol, and beads are dispensed into specified locations on a temperature module and magnetic rack, preparing for subsequent steps in DNA purification and processing.

from opentrons.types import Point
metadata = {
    'protocolName': 'mtb',
    'author': 'MGI-X',
    'apiLevel': '2.14'
}

def run(ctx):
    #Reagent loading
    temp_rack["A1"].load_liquid(liquid={round}nd_TE_buffer, volume = 2nd_purification_TE_Buffer_tolerance)
    mag_rack["A1"].load_liquid(liquid={round}nd_pcr_product, volume = 25)
    plate_1["A1"].load_liquid(liquid=Ethonal, volume = 480)
    temp_rack["B1"].load_liquid(liquid=Beads, volume = 50)