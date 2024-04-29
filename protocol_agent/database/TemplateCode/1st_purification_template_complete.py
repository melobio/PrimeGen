#ALL Code description: This script is suitable for the first round of purification steps in a two-step experiment. This script configures the Opentrons OT-2 robot for a DNA purification protocol using magnetic beads. It initializes hardware modules, defines reagents, and executes a series of pipetting tasks for mixing, binding, washing, and eluting DNA.**Module Initialization:** It loads a temperature module, a magnetic module, and appropriate labware for reagents and samples.**Reagent Definition:** Liquids such as TE buffer, ethanol, magnetic beads, and PCR products are defined with specific properties like color and descriptions.**Reagent Loading:** Each reagent is loaded into specified wells on the labware for use in the protocol.**Pipetting and Purification:** The script performs precise pipetting actions to mix magnetic beads with the PCR product, wash bound DNA with ethanol, and finally elute the purified DNA in TE buffer.**Magnetic Operations:** The magnetic module is engaged to separate the beads from the supernatant during cleaning steps, and disengaged after the final elution step.This automated process enhances the efficiency and reproducibility of DNA purification tasks in molecular biology experiments.
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
    {round}st_pcr_product = ctx.define_liquid(
        name=f"{round}st_PCR_product",
        description=f"{round}st_pcr_product",
        display_color="#00FF00",
        )
    TE_buffer= ctx.define_liquid(
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
    temp_rack["A1"].load_liquid(liquid=TE_buffer, volume = 13)
    mag_rack["A1"].load_liquid(liquid={round}st_pcr_product, volume = 25)
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
        volume=1st_purification_TE_Buffer_elute_DNA,
        source=temp_rack.wells_by_name()["A1"],
        dest=mag_rack.wells_by_name()["A1"],
        new_tip='once',
        )  
    
    #Remove magnetism
    mag_mod.disengage()
