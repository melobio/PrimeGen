#Code block description: This script block is suitable for the first round of PCR reaction steps in a two-step experiment. This script sets up a PCR (Polymerase Chain Reaction) protocol on the Opentrons OT-2 robot, specifically configuring the thermocycler module (tc_mod). It involves closing the lid, setting the lid temperature, and managing various temperature stages critical for PCR. The script controls the block temperature for initial denaturation, annealing, and extension phases, and it executes a detailed temperature profile multiple times according to specified PCR cycles. After the cycling is complete, the script cools down the system to end the reaction, ensuring that samples are maintained at a low temperature post-PCR.

from opentrons.types import Point
metadata = {
    'protocolName': 'mtb',
    'author': 'MGI-X',
    'apiLevel': '2.14'
}

def run(ctx):
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
