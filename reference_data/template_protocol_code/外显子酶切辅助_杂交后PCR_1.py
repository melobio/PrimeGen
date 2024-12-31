from opentrons import protocol_api

metadata = {
    'protocolName': '外显子组酶切辅助_杂交后PCR_1',
    'author': 'MGIX',
    'description': 'Protocol for DNA fragmentation and recovery',
    'apiLevel': '2.14'
}

# Transfer function
def transfer_all(pipette, source_plate, dest_plate, transfer_info, mix_flag=0, mix_value=[]):
    """
    General function for transferring samples.
    Args:
    - pipette: The pipette to use
    - source_plate: Source plate
    - dest_plate: Destination plate
    - transfer_info: List containing transfer information, each element is a tuple (source column name, destination column name, volume)
    - mix_value: Mixing information, a list [times, volume]
    """
    source_map = source_plate.columns_by_name()
    dest_map = dest_plate.columns_by_name()
    for source_col, dest_col, volume in transfer_info:
        temp_source = source_map.get(source_col)
        temp_dest = dest_map.get(dest_col)
        pipette.pick_up_tip()
        pipette.transfer(volume, temp_source, temp_dest, new_tip='never')
        if mix_flag == 1:
            for well in temp_dest:
                pipette.mix(mix_value[0], mix_value[1], well)
        pipette.drop_tip()

def run(protocol: protocol_api.ProtocolContext):
    # Load 200ul tip rack
    tips200_1 = protocol.load_labware('opentrons_96_filtertiprack_200ul', '2')
    temp_module = protocol.load_module('temperature module', '3')
    temp_plate = temp_module.load_labware('biorad_96_wellplate_200ul_pcr')
    # Load thermocycler module, ensure to specify a correct slot number
    thermocycler_module = protocol.load_module('thermocycler', '7')
    thermocycler_plate = thermocycler_module.load_labware('biorad_96_wellplate_200ul_pcr')
    # Pipettes
    right_pipette = protocol.load_instrument('p300_multi_gen2', mount='right', tip_racks=[tips200_1])
    # Define liquids
    # labeling liquids in wells
    sample = protocol.define_liquid(
        name="sample",
        description="sample",
        display_color="#00FF00",
    )
    PCR_MIX = protocol.define_liquid(
        name="PCR_MIX",
        description="PCR_MIX",
        display_color="#00FF00",
    )
    # Load liquids
    for well in temp_plate.columns_by_name()['1']:
        temp_plate.wells_by_name()[well.display_name].load_liquid(liquid=PCR_MIX, volume=182)
    for well in temp_plate.columns_by_name()['2']:
        temp_plate.wells_by_name()[well.display_name].load_liquid(liquid=PCR_MIX, volume=182)
    for well in temp_plate.columns_by_name()['3']:
        temp_plate.wells_by_name()[well.display_name].load_liquid(liquid=PCR_MIX, volume=182)
    for well in temp_plate.columns_by_name()['4']:
        temp_plate.wells_by_name()[well.display_name].load_liquid(liquid=PCR_MIX, volume=182)
    well_name_list = [well.display_name for well in thermocycler_plate.wells()]
    for well_name in well_name_list:
        thermocycler_plate.wells_by_name()[well_name].load_liquid(liquid=sample, volume=19)

    # Execute instructions
    # 1. Set the temperature module to 4 degrees and pause to wait for the set temperature
    temp_module.set_temperature(4)
    temp_module.await_temperature(4)
    # Open thermocycler
    thermocycler_module.open_lid()
    thermocycler_module.set_block_temperature(4)
    # Transfer
    transfer_times = 12
    transfer_info_right1 = []
    for ii in range(transfer_times):
        if ii < 3:
            transfer_info_right1.append(('1', str(ii + 1), 56))
        elif ii >= 3 and ii < 6:
            transfer_info_right1.append(('2', str(ii + 1), 56))
        elif ii >= 6 and ii < 9:
            transfer_info_right1.append(('3', str(ii + 1), 56))
        else:
            transfer_info_right1.append(('4', str(ii + 1), 56))
    transfer_all(right_pipette, temp_plate, thermocycler_plate, transfer_info_right1, mix_flag=1, mix_value=[5, 25])

    thermocycler_module.close_lid()
    thermocycler_module.set_lid_temperature(105)
    steps1 = [{'temperature': 95, 'hold_time_minutes': 3}]
    steps2 = [{'temperature': 98, 'hold_time_seconds': 20}, {'temperature': 60, 'hold_time_seconds': 15},
              {'temperature': 72, 'hold_time_seconds': 30}]
    steps3 = [{'temperature': 72, 'hold_time_minutes': 10}, {'temperature': 4, 'hold_time_minutes': 1}]
    thermocycler_module.execute_profile(steps=steps1, repetitions=1)
    thermocycler_module.execute_profile(steps=steps2, repetitions=12)
    thermocycler_module.execute_profile(steps=steps3, repetitions=1)
    thermocycler_module.set_block_temperature(4, hold_time_seconds=None)
    thermocycler_module.deactivate_lid()
    thermocycler_module.open_lid()
