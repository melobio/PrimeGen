#ALL Code description: This script is suitable for the make DNB steps in a two-step experiment.This Opentrons OT-2 script automates the setup for a complex DNA library preparation protocol involving the creation of DNA NanoBalls (DNBs). The script initializes temperature control and thermocycler modules, loads reagents such as Low TE buffer, dsDNA library, DNB buffers, and enzyme mixes onto a temperature module. Detailed pipetting steps transfer these reagents into a thermocycler for precise thermal cycling during DNB formation and rolling circle amplification stages. This includes temperature steps for primer hybridization and amplification, followed by adding a stop buffer to terminate the reaction, ensuring the process is accurate and reproducible for next-generation sequencing preparation.

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
    # pipette and tiprack
    tiprack = [ctx.load_labware('opentrons_96_filtertiprack_20ul', slot) for slot in ['3']]
    pip20_single = ctx.load_instrument('p20_single_gen2', 'right', tip_racks=[*tiprack])
    pip20_multi = ctx.load_instrument('p20_multi_gen2', 'left', tip_racks=[*tiprack])

    #Reagent definition
    Low_TE_buffer = ctx.define_liquid(
        name="Low_TE_buffer",
        description="Low_TE_buffer",#Low_TE_buffer
        display_color="#00FF00",
        )
    dsDNA_libarary = ctx.define_liquid(
        name="dsDNA_libarary",
        description="dsDNA_libarary",#dsDNA_libarary
        display_color="#00FF00",
        )
    make_DNB_buffer = ctx.define_liquid(
        name="make_DNB_buffer",
        description="make_DNB_buffer",#make DNB buffer
        display_color="#00FF00",
        )
    make_DNB_enzyme_mix_1= ctx.define_liquid(
        name="make_DNB_enzyme_mix_1",
        description="make_DNB_enzyme_mix_1",#make_DNB_enzyme_mix_1
        display_color="#00FF00",
        )
    make_DNB_enzyme_mix_2= ctx.define_liquid(
        name="make_DNB_enzyme_mix_2",
        description="make_DNB_enzyme_mix_2",#make_DNB_enzyme_mix_2
        display_color="#00FF00",
        )
    make_DNB_stop_buffer = ctx.define_liquid(
        name="Beads",
        description="Beads",#Terminate DNB reaction
        display_color="#00FF00",
        )
    
    #Reagent loading
    #Users need to upload Qubit’s concentration information for further calculations
    temp_rack["A1"].load_liquid(liquid=Low_TE_buffer, volume = Low_TE_buffer_Volume_tolerance)
    temp_rack["B1"].load_liquid(liquid=dsDNA_libarary, volume = dsDNA_library_volume_tolerance)
    temp_rack["C1"].load_liquid(liquid=make_DNB_buffer, volume = Make_DNB_Buffer_Volume_tolerance)
    temp_rack["D1"].load_liquid(liquid=make_DNB_enzyme_mix_1, volume = Make_DNB_Enzyme_Mix_I_tolerance)
    temp_rack["E1"].load_liquid(liquid=make_DNB_enzyme_mix_2, volume = Make_DNB_Enzyme_Mix_II_tolerance)
    temp_rack["F1"].load_liquid(liquid=make_DNB_stop_buffer, volume = Make_DNB_Stop_Buffer_tolerance)
    
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

    #Make DNB reaction 1
    tc_mod.close_lid()
    tc_mod.set_lid_temperature(temperature=105)
    tc_mod.set_block_temperature(temperature=Primer_hybridization_Reaction_step1_temperature,hold_time_minutes=3,block_max_volume=40)
    tc_mod.set_block_temperature(temperature=Primer_hybridization_Reaction_step2_temperature,hold_time_minutes=3,block_max_volume=40)
    tc_mod.set_block_temperature(temperature=Primer_hybridization_Reaction_step3_temperature,hold_time_minutes=3,block_max_volume=40)

    #Pipetting
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
    
    #PCR reaction configuration
    tc_mod.close_lid()
    tc_mod.set_lid_temperature(temperature=37)
    tc_mod.set_block_temperature(temperature=Rolling circle amplification_step1_temperature,hold_time_minutes=25,block_max_volume=85)
    tc_mod.set_block_temperature(temperature=Rolling circle amplification_step2_temperature,hold_time_minutes=3,block_max_volume=85)
    tc_mod.open_lid()
    
    #Pipetting
    pip20_single.flow_rate.aspirate = 3.78
    pip20_single.flow_rate.dispense = 3.78
    for well_name in ['A1']:
        pip20_single.transfer(
            volume=Make_DNB_Stop_Buffer,
            source=temp_rack.wells_by_name()["F1"],
            dest=tc_rack.wells_by_name()[well_name],
            new_tip='always',
            ) 