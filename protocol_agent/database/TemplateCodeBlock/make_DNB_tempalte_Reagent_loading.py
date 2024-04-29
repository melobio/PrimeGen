#Code block description: This script block is suitable for the make DNB steps in a two-step experiment. This script automates the loading of specific reagents into a temperature-controlled rack on the Opentrons OT-2 robot. It includes Low TE buffer, dsDNA library, DNB buffer, two types of DNB enzyme mixes, and a stop buffer, each loaded according to specified tolerances necessary for DNA library preparation.

from opentrons.types import Point
metadata = {
    'protocolName': 'mtb',
    'author': 'MGI-X',
    'apiLevel': '2.14'
}

def run(ctx):
    #Reagent loading
    #Users need to upload Qubit’s concentration information for further calculations
    temp_rack["A1"].load_liquid(liquid=Low_TE_buffer, volume = Low_TE_buffer_Volume_tolerance)
    temp_rack["B1"].load_liquid(liquid=dsDNA_libarary, volume = dsDNA_library_volume_tolerance)
    temp_rack["C1"].load_liquid(liquid=make_DNB_buffer, volume = Make_DNB_Buffer_Volume_tolerance)
    temp_rack["D1"].load_liquid(liquid=make_DNB_enzyme_mix_1, volume = Make_DNB_Enzyme_Mix_I_tolerance)
    temp_rack["E1"].load_liquid(liquid=make_DNB_enzyme_mix_2, volume = Make_DNB_Enzyme_Mix_II_tolerance)
    temp_rack["F1"].load_liquid(liquid=make_DNB_stop_buffer, volume = Make_DNB_Stop_Buffer_tolerance)