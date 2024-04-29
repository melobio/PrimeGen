#Code block description: This script block is suitable for the second round of PCR reaction steps in a two-step experiment. This script configures a PCR protocol on the Opentrons OT-2 using a thermocycler module. It sets specific temperatures for the lid and block, runs a detailed thermal profile for DNA amplification, and maintains low temperatures post-PCR to stabilize the products.

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
    tc_mod.set_block_temperature(temperature=37,hold_time_minutes=5,block_max_volume=2nd_tcmod_block_max_volume)
    tc_mod.set_block_temperature(temperature=95,hold_time_minutes=10,block_max_volume=2nd_tcmod_block_max_volume)
    profile = [
        {'temperature':95, 'hold_time_seconds':20},
        {'temperature':64, 'hold_time_seconds':30},
        {'temperature':60, 'hold_time_seconds':30},
        {'temperature':72, 'hold_time_seconds':30},
        ]
    tc_mod.execute_profile(steps=profile, repetitions=2nd_PCR_Cycles, block_max_volume=2nd_tcmod_block_max_volume)
    tc_mod.set_block_temperature(temperature=4,hold_time_minutes=3,block_max_volume=2nd_tcmod_block_max_volume) 
    tc_mod.set_lid_temperature(temperature=105)
    tc_mod.set_block_temperature(temperature=4) 
