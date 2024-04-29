#Code block description: This script block is suitable for the make DNB steps in a two-step experiment. This script manages the pipetting tasks for a DNA library preparation on the Opentrons OT-2. It sets precise temperatures, opens the thermocycler lid, and adjusts pipetting parameters. Low TE buffer, dsDNA library, and DNB buffer are meticulously transferred and mixed within specified wells, ensuring precise and controlled sample preparation.

from opentrons.types import Point
metadata = {
    'protocolName': 'mtb',
    'author': 'MGI-X',
    'apiLevel': '2.14'
}

def run(ctx):
    #Pipetting
    temp_mod.set_temperature(celsius=4) 
    tc_mod.open_lid()
    tc_mod.set_block_temperature(temperature=4)
    pip20_single.well_bottom_clearance.aspirate = 1  
    pip20_single.well_bottom_clearance.dispense = 0.5  
    pip20_single.flow_rate.aspirate = 3.78
    pip20_single.flow_rate.dispense = 20
    
    #Pipetting Low_TE_buffer
    for well_name in ['A1']:
        pip20_single.transfer(
            volume=Low_TE_buffer_Volume,
            source=temp_rack.wells_by_name()["A1"],
            dest=tc_rack.wells_by_name()[well_name],
            new_tip='always',
            )
    #Pipetting dsDNA_library
    for well_name in ['A1']:
        pip20_single.transfer(
            volume=dsDNA_library_volume,
            source=temp_rack.wells_by_name()["B1"],
            dest=tc_rack.wells_by_name()[well_name],
            new_tip='always',
            )
    #Pipetting Make_DNB_Buffer and mix
    for well_name in ['A1']:
        pip20_single.transfer(
            volume=Make_DNB_Buffer_Volume,
            source=temp_rack.wells_by_name()["C1"],
            dest=tc_rack.wells_by_name()[well_name],
            new_tip='always',
            mix_after=(5, 15),
            )  