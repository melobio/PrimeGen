#ALL Code description: This script is suitable for the first round of PCR reaction steps in a two-step experiment. This comprehensive Opentrons OT-2 script outlines a complete workflow for a PCR experiment. It includes module loading, reagent definition, reagent loading, pipetting operations, and PCR reaction configurations. The script manages different modules like temperature modules and thermocyclers, defines and loads various reagents such as DNA templates, PCR enzyme mixes, and primer pools, and handles liquid transfers with precision. Pipetting parameters are set for optimal liquid handling. Finally, a detailed PCR cycling profile is executed, ensuring the accurate thermal cycling required for successful PCR amplification.

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
    # Plates
    plate_1 = ctx.load_labware('biorad_96_wellplate_200ul_pcr', '1', 'reg1') 
    # The number of tipracks needs to be calculated
    tiprack = [ctx.load_labware('opentrons_96_filtertiprack_20ul', slot) for slot in 1st_Tip_rack_number] 
    pip20_single = ctx.load_instrument('p20_single_gen2', 'right', tip_racks=[*tiprack]) 
    pip20_multi = ctx.load_instrument('p20_multi_gen2', 'left', tip_racks=[*tiprack]) 

 
    #Reagent definition
    DNA_TEMPLATE = ctx.define_liquid(
    name=f"DNA_TEMPLATE_{PCR_round}st_PCR",
    description=f"DNA_TEMPLATE_{PCR_round}st_PCR", 
    display_color="#00FF00",
    )
    Barcode_{PCR_round} = ctx.define_liquid(
    name=f"Barcode_{PCR_round}",
    description="Barcode_{PCR_round}", 
    display_color="#00FF00",
    )
    primer_pools = ctx.define_liquid(
    name="primer pools",
    description="primer pools", 
    display_color="#00FF00", 
    )
    PCR_Enzyme_Mix = ctx.define_liquid(
    name=f"{PCR_round}st_PCR_Enzyme_Mix",
    description=f"{PCR_round}st PCR Enzyme Mix", 
    display_color="#00FF00",
    )

    #Reagent loading
    for well in plate_1.wells()[:96]:
        well.load_liquid(liquid=Barcode_{PCR_round}, volume=30) 

    for well_name in {1st_PCR_Mix_well}：
         temp_rack[well_name].load_liquid(liquid={PCR_round}st_PCR_Enzyme_Mix, volume=1st_PCR_MIX_Volume_fault_tolerance_volume) 
    
    temp_rack[temp_rack.wells()[len(1st_PCR_Mix_well)+1].load_liquid(liquid=primer_pools, volume=Primer_Pool_well_Volume) 
    
    for well_name in DNA_TEMPLATE_WELL: 
        tc_rack.wells_by_name()[well_name].load_liquid(liquid=DNA_TEMPLATE, volume=2.5) 

    #Pipetting
    temp_mod.set_temperature(celsius=4) #The temperature module cools down to 4℃
    tc_mod.open_lid() #thermocyclerModule opening lid
    tc_mod.set_block_temperature(temperature=4) #thermocyclerModule block(reaction well)cools down to 4℃
    
    pip20_single.flow_rate.aspirate = 3.78 #Aspiration speed
    pip20_single.flow_rate.dispense = 20 #dispense speed
    #PCR MIX Pipetting
    pcr_mix_wells = temp_rack.wells()[:len(1st_PCR_Mix_well)]  
    for well in tc_rack.wells()[: sample size]:  
        source_well = pcr_mix_wells[int(well.get_name()[1:]) //(96_well_volume//1st_PCR_MIX_Volume_fault_tolerance_volume))]  
        pip20.transfer(volume=1st_PCR_MIX_Volume, 
            source=source_well, 
            dest=well, 
            new_tip='always')  
    #DNA TEMPLATE Pipetting
    pip20_single.well_bottom_clearance.aspirate = 1  
    pip20_single.well_bottom_clearance.dispense = 0.5  
    source_plate = [plate_1[i] for i in  DNA_TEMPLATE_WELL]  
    dest_plate = [tc_rack[i] for i in  DNA_TEMPLATE_WELL]  
    pip20_single.transfer(    
            volume=2,    
            source=source_plate ,    
            dest=dest_plate ,    
            new_tip='always',    
            )    
    #Primer Pool Pipetting
    for well_name in DNA_TEMPLATE_WELL:
        pip20_single.pick_up_tip()
        pip20_single.aspirate(Primer_Pool_Volume, temp_rack.wells_by_name()[temp_rack.wells()[len(1st_PCR_Mix_well)+1]].bottom(z=1))
        pip20_single.dispense(Primer_Pool_Volume, tc_rack.wells_by_name()[well_name].bottom(z=0.5))
        ctx.delay(seconds=1)
        pip20_single.drop_tip()
    ctx.pause('Seal the Bio-Rad 96 Well Plate and centrifuge, then put it back into the thermocycler Module and continue the experiment.')

    #PCR reaction configuration
    tc_mod.close_lid()
    tc_mod.set_lid_temperature(temperature=105)
    tc_mod.set_block_temperature(temperature=37,hold_time_minutes=5,block_max_volume=1st_tcmod_block_max_volume)
    tc_mod.set_block_temperature(temperature=95,hold_time_minutes=10,block_max_volume=1st_tcmod_block_max_volume)
    profile = [
        {'temperature':95, 'hold_time_seconds':20},
        {'temperature':64, 'hold_time_seconds':30},
        {'temperature':60, 'hold_time_seconds':30},
        {'temperature':72, 'hold_time_seconds':30},
        ]
    tc_mod.execute_profile(steps=profile, repetitions=1st_PCR_Cycles, block_max_volume=1st_tcmod_block_max_volume)
    tc_mod.set_block_temperature(temperature=4,hold_time_minutes=3,block_max_volume=1st_tcmod_block_max_volume) 
    tc_mod.set_lid_temperature(temperature=105)
    tc_mod.set_block_temperature(temperature=4) 
