#Code block description: This script block is suitable for the first round of purification steps in a two-step experiment. This script automates DNA purification using magnetic beads on an Opentrons OT-2 robot. It includes precise pipetting for mixing, washing with ethanol, and DNA elution in TE buffer. The process involves magnetic separation, timed delays for settling, and careful liquid transfer to ensure efficient purification.
from opentrons.types import Point
metadata = {
    'protocolName': 'mtb',
    'author': 'MGI-X',
    'apiLevel': '2.14'
}

def run(ctx):
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
