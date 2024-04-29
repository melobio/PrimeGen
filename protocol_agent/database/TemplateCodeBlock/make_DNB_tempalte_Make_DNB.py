#Code block description: This script block is suitable for the make DNB steps in a two-step experiment. This script controls the Opentrons OT-2 for a DNA NanoBall (DNB) creation process, involving precise thermal cycling for primer hybridization and rolling circle amplification. It handles specific pipetting sequences for enzyme mixes and stop buffer into a thermocycler module, ensuring accurate reagent delivery and mixing for effective DNA amplification.

from opentrons.types import Point
metadata = {
    'protocolName': 'mtb',
    'author': 'MGI-X',
    'apiLevel': '2.14'
}

def run(ctx):
    #Make DNB
    ##Make DNB reaction 1
    tc_mod.close_lid()
    tc_mod.set_lid_temperature(temperature=105)
    tc_mod.set_block_temperature(temperature=Primer_hybridization_Reaction_step1_temperature,hold_time_minutes=3,block_max_volume=40)
    tc_mod.set_block_temperature(temperature=Primer_hybridization_Reaction_step2_temperature,hold_time_minutes=3,block_max_volume=40)
    tc_mod.set_block_temperature(temperature=Primer_hybridization_Reaction_step3_temperature,hold_time_minutes=3,block_max_volume=40)

    ##Pipetting
    tc_mod.open_lid()
    for well_name in ['A1']:
        pip20_single.transfer(
            volume=Make_DNB_Enzyme_Mix_I,
            source=temp_rack.wells_by_name()["D1"],
            dest=tc_rack.wells_by_name()[well_name],
            new_tip='always',
            mix_after=(5, 15),
            )  
        
    for well_name in ['A1']:
        pip20_single.transfer(
            volume=Make_DNB_Enzyme_Mix_II,
            source=temp_rack.wells_by_name()["E1"],
            dest=tc_rack.wells_by_name()[well_name],
            new_tip='always',
            mix_after=(5, 15),
            )  
    
    ##PCR reaction configuration
    tc_mod.close_lid()
    tc_mod.set_lid_temperature(temperature=37)
    tc_mod.set_block_temperature(temperature=Rolling circle amplification_step1_temperature,hold_time_minutes=25,block_max_volume=85)
    tc_mod.set_block_temperature(temperature=Rolling circle amplification_step2_temperature,hold_time_minutes=3,block_max_volume=85)
    tc_mod.open_lid()
    
    ##Pipetting
    pip20_single.flow_rate.aspirate = 3.78
    pip20_single.flow_rate.dispense = 3.78
    for well_name in ['A1']:
        pip20_single.transfer(
            volume=Make_DNB_Stop_Buffer,
            source=temp_rack.wells_by_name()["F1"],
            dest=tc_rack.wells_by_name()[well_name],
            new_tip='always',
            ) 