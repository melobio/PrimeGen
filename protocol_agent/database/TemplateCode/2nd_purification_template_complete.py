#ALL Code description: This script is suitable for the second round of purification steps in a two-step experiment.This script orchestrates a DNA purification process on the Opentrons OT-2 robot, utilizing temperature, magnetic modules, and pipetting instruments. It involves defining and loading reagents like PCR products, TE buffer, ethanol, and magnetic beads into specific labware. Key operations include transferring and mixing magnetic beads with PCR products for DNA binding, followed by sequential washing steps with ethanol to clean the beads. After drying the beads, TE buffer is added to elute the DNA, which is then mixed to dissolve properly. The script concludes with magnetic separation to collect the purified DNA, ensuring thorough, automated purification suitable for downstream applications.
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
    mag_mod = ctx.load_module('magnetic module gen2', '3')
    mag_rack = mag_mod.load_labware('biorad_96_wellplate_200ul_pcr')
    # Plates
    plate_1 = ctx.load_labware('nest_96_wellplate_2ml_deep', '1', 'reg1')
    # pipette and tiprack
    tiprack = [ctx.load_labware('opentrons_96_filtertiprack_20ul', slot) for slot in ['2']] 
    pip20_single = ctx.load_instrument('p20_single_gen2', 'right', tip_racks=[*tiprack]) 
    pip20_multi = ctx.load_instrument('p20_multi_gen2', 'left', tip_racks=[*tiprack]) 


    #Reagent definition
    {round}nd_pcr_product = ctx.define_liquid(
        name=f"{round}nd_PCR_product",
        description=f"{round}nd_pcr_product",
        display_color="#00FF00",
        )
    {round}nd_TE_buffer= ctx.define_liquid(
        name="TE buffer",
        description="TE buffer",#Used to dissolve DNA
        display_color="#00FF00",
        )
    Ethonal = ctx.define_liquid(
        name="Ethonal",
        description="Ethonal",#Used to clean magnetic beads
        display_color="#00FF00",
        )
    Beads = ctx.define_liquid(
        name="Beads",
        description="Beads",#Used to bind DNA
        display_color="#00FF00",
        )

    #Reagent loading
    temp_rack["A1"].load_liquid(liquid={round}nd_TE_buffer, volume = 2nd_purification_TE_Buffer_tolerance)
    mag_rack["A1"].load_liquid(liquid={round}nd_pcr_product, volume = 25)
    plate_1["A1"].load_liquid(liquid=Ethonal, volume = 480)
    temp_rack["B1"].load_liquid(liquid=Beads, volume = 50)

    #Pipetting and purification
    pip20_single.well_bottom_clearance.dispense = 0.5  
    pip20_single.transfer(
        volume=25,
        source=temp_rack.wells_by_name()["B1"],
        dest=mag_rack.wells_by_name()["A1"],
        new_tip='once',
        )
    
    #Magnetic beads mixed
    pip20_single.well_bottom_clearance.aspirate = 0.5  
    pip20_single.pick_up_tip()
    pip20_single.mix(repetitions=10, volume=15,mag_rack['A3'])
    pip20_single.mix(10, 15, mag_rack['A1'])    
    pip20_single.drop_tip()  
    ctx.delay(minutes=5)
    mag_mod.engage(height_from_base=5)
    ctx.delay(minutes=5)
    
    #Remove supernatant
    pip20_single.well_bottom_clearance.aspirate = 1  
    pip20_single.well_bottom_clearance.dispense = 0.5  
    pip20_single.transfer(
        volume=50,
        source=mag_rack.wells_by_name()["A1"],
        dest=temp_rack.wells_by_name()["C1"], #waste hole
        new_tip='once',
        )
    
    #Ethanol cleaning magnetic beads
    pip20_single.transfer(
        volume=160,
        source=plate_1.wells_by_name()["A1"],
        dest=mag_rack.wells_by_name()["A1"],
        new_tip='once',
        )
    ctx.delay(seconds=30)
    
    #Remove Ethanol 
    pip20_single.transfer(
        volume=160,
        source=mag_rack.wells_by_name()["A1"],
        dest=temp_rack.wells_by_name()["D1"], #waste hole
        new_tip='once',
        )  
    
    #Repeat cleaning
    pip20_single.transfer(
        volume=160,
        source=plate_1.wells_by_name()["A1"],
        dest=mag_rack.wells_by_name()["A1"],
        new_tip='once',
        )
    ctx.delay(seconds=30)
    pip20_single.transfer(
        volume=160,
        source=mag_rack.wells_by_name()["A1"],
        dest=temp_rack.wells_by_name()["E1"],
        new_tip='once',
        )  
    
    #dry ethanol
    ctx.delay(minutes=5)
    
    #TE Buffer dissolves DNA
    pip20_single.transfer(
        volume=2nd_purification_TE_Buffer_elute_DNA,
        source=temp_rack.wells_by_name()["A1"],
        dest=mag_rack.wells_by_name()["A1"],
        new_tip='once',
        )  
    
    #Remove magnetism
    mag_mod.disengage()

    #Mix well to dissolve the DNA
    pip20_single.well_bottom_clearance.aspirate = 0.5  
    pip20_single.pick_up_tip()
    pip20_single.mix(10, 15, mag_rack['A1'])    
    pip20_single.drop_tip()  

    ctx.delay(minutes=5)

    mag_mod.engage(height_from_base=5)    
    ctx.delay(minutes=5)
    
    pip20_single.transfer(
        volume=2nd_purification_supernatant_volumen,
        source=mag_rack.wells_by_name()["A1"],
        dest=temp_rack.wells_by_name()["A2"],
        new_tip='once',
        )  