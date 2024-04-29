#Code block description: This script block is suitable for the second round of PCR reaction steps in a two-step experiment. This script configures pipetting operations for a PCR protocol on the Opentrons OT-2. It sets module temperatures, manages pipette flow rates and clearances, and precisely transfers PCR mixes, DNA templates, and other reagents between various wells. The script includes pauses for manual interventions like sealing and centrifugation.


from opentrons.types import Point
metadata = {
    'protocolName': 'mtb',
    'author': 'MGI-X',
    'apiLevel': '2.14'
}

def run(ctx):
    #Pipetting
    temp_mod.set_temperature(celsius=4) #The temperature module cools down to 4℃
    tc_mod.open_lid() #thermocyclerModule opening lid
    tc_mod.set_block_temperature(temperature=4) #thermocyclerModule block(reaction well)cools down to 4℃
    pip20_single.flow_rate.aspirate = 3.78 #Aspiration speed
    pip20_single.flow_rate.dispense = 20 #dispense speed
    pcr_mix_wells = temp_rack.wells()[:1]  
    for well in tc_rack.wells()[:1]:  
        source_well = source_well = pcr_mix_wells[int(well.get_name()[1:]) //(96_well_volume//2nd_PCR_MIX_Volume_fault_tolerance_volume))]  
        pip20.transfer(2nd_PCR_MIX_Volume, source_well, well, new_tip='always')  
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
    for well_name in DNA_TEMPLATE_WELL:
        pip20_single.pick_up_tip()
        pip20_single.aspirate(U-R_well_Volume, temp_rack.wells_by_name()[temp_rack.wells()[2]].bottom(z=1))
        pip20_single.dispense(U-R_well_Volume, tc_rack.wells_by_name()[well_name].bottom(z=0.5))
        ctx.delay(seconds=1)
        pip20_single.drop_tip()
    ctx.pause('Seal the Bio-Rad 96 Well Plate and centrifuge, then put it back into the thermocycler Module and continue the experiment.')
