from opentrons import protocol_api

metadata = {
    'protocolName': '单细胞染色体拷贝数变异_酶切打断_2',
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
        if mix_flag == 1:
            pipette.transfer(volume, temp_source, temp_dest, new_tip='never', mix_after=(mix_value[0], mix_value[1]))
        else:
            pipette.transfer(volume, temp_source, temp_dest, new_tip='never')
        pipette.drop_tip()


def calculate_liquid(sample_num, liquid_volume):
    liquid_list = []
    if sample_num < 8:
        per_sample_liquid = liquid_volume / sample_num
        temp_num = sample_num
        for ii in range(8):
            if temp_num > 0:
                liquid_list.append(per_sample_liquid)
                temp_num = temp_num - 1
            else:
                liquid_list.append(0.0)
        return liquid_list
    else:
        count1 = int(sample_num / 8)
        count2 = sample_num % 8
        per_sample_liquid = liquid_volume / (8 * (count1 + 1) + count2)

        for ii in range(8):
            if count2 > 0:
                liquid_list.append((count1 + 2) * per_sample_liquid)
                count2 = count2 - 1
            else:
                liquid_list.append((count1 + 1) * per_sample_liquid)
        return liquid_list


def run(protocol: protocol_api.ProtocolContext):
    Sample_nums = input_sample_nums_1
    Fragment_MIX_volume = input_Fragment_MIX_volume_2
    TE_Buffer_volume = input_TE_Buffer_volume_2

    use_col_nums = int(Sample_nums / 8)
    if Sample_nums % 8 != 0:
        use_col_nums = use_col_nums + 1

    # Load 20ul tip rack
    tips20_1 = protocol.load_labware('opentrons_96_tiprack_20ul', '4')
    tips20_2 = protocol.load_labware('opentrons_96_tiprack_20ul', '5')
    # Load 200ul filter tip rack
    tips200_1 = protocol.load_labware('opentrons_96_tiprack_300ul', '2')
    plate1 = protocol.load_labware('biorad_96_wellplate_200ul_pcr', '1')
    temp_module = protocol.load_module('temperature module', '3')
    temp_plate = temp_module.load_labware('biorad_96_wellplate_200ul_pcr')
    # Load thermocycler module, ensure to specify a correct slot number
    thermocycler_module = protocol.load_module('thermocycler', '7')
    thermocycler_plate = thermocycler_module.load_labware('biorad_96_wellplate_200ul_pcr')
    # Pipettes
    left_pipette = protocol.load_instrument('p20_multi_gen2', mount='left', tip_racks=[tips20_1,tips20_2])
    right_pipette = protocol.load_instrument('p300_multi_gen2', mount='right', tip_racks=[tips200_1])
    left_pipette.flow_rate.aspirate = 3.78  # Set aspirate speed to 50 μL/s
    left_pipette.flow_rate.dispense = 7.56
    right_pipette.flow_rate.aspirate = 46.43  # Set aspirate speed to 50 μL/s
    right_pipette.flow_rate.dispense = 92.86
    # Define liquids
    # labeling liquids in wells
    sample = protocol.define_liquid(
        name="sample",
        description="sample",
        display_color="#00FF00",
    )
    Fragment_MIX = protocol.define_liquid(
        name="Fragment_MIX",
        description="Fragment_MIX",
        display_color="#00FF00",
    )
    TE_Buffer = protocol.define_liquid(
        name="TE_Buffer",
        description="TE_Buffer",
        display_color="#00FF00",
    )
    # Load liquids
    all_plate_well = [f"{chr(65 + i)}{j}" for j in range(1, 12 + 1) for i in range(8)]
    liquid_list = calculate_liquid(Sample_nums, Fragment_MIX_volume)
    for ii in range(8):
        well_name = all_plate_well[ii]
        temp_plate.wells_by_name()[well_name].load_liquid(liquid=Fragment_MIX, volume=liquid_list[ii])

    if use_col_nums > 6:
        for ii in range(8):
            temp_plate.wells_by_name()[all_plate_well[8:16][ii]].load_liquid(liquid=TE_Buffer, volume=130)
        liquid_list2 = calculate_liquid(Sample_nums - 48, TE_Buffer_volume - 130)
        for ii in range(8):
            temp_plate.wells_by_name()[all_plate_well[16:24][ii]].load_liquid(liquid=TE_Buffer, volume=liquid_list2[ii])
    else:
        liquid_list3 = calculate_liquid(Sample_nums, TE_Buffer_volume)
        for ii in range(8):
            temp_plate.wells_by_name()[all_plate_well[8:16][ii]].load_liquid(liquid=TE_Buffer, volume=liquid_list3[ii])
    use_plate_well = all_plate_well[:Sample_nums]
    for well_name in use_plate_well:
        plate1.wells_by_name()[well_name].load_liquid(liquid=sample, volume=60)

    # Execute instructions
    # 1. Set the temperature module to 4 degrees and pause to wait for the set temperature
    temp_module.set_temperature(4)
    temp_module.await_temperature(4)
    # Start thermocycler
    thermocycler_module.open_lid()
    thermocycler_module.set_block_temperature(4)
    # Transfer
    transfer_times = use_col_nums
    transfer_info_right1 = []
    for ii in range(transfer_times):
        transfer_info_right1.append((str(ii + 1), str(ii + 1), 24))

    transfer_all(right_pipette, plate1, thermocycler_plate, transfer_info_right1)

    transfer_info_left2 = []
    for ii in range(transfer_times):
        transfer_info_left2.append(('1', str(ii + 1), 6))

    transfer_all(left_pipette, temp_plate, thermocycler_plate, transfer_info_left2, mix_flag=1, mix_value=[5, 10])

    thermocycler_module.close_lid()
    thermocycler_module.set_lid_temperature(105)
    steps1 = [{'temperature': 37, 'hold_time_seconds': 300}, {'temperature': 75, 'hold_time_seconds': 900},
              {'temperature': 4, 'hold_time_seconds': 60}]
    thermocycler_module.execute_profile(steps=steps1, repetitions=1, block_max_volume=60)
    thermocycler_module.set_block_temperature(4, hold_time_seconds=None)
    thermocycler_module.deactivate_lid()
    thermocycler_module.open_lid()

    transfer_info_left3 = []
    for ii in range(transfer_times):
        if ii<6:
            transfer_info_left3.append(('2', str(ii + 1), 20))
        else:
            transfer_info_left3.append(('3', str(ii + 1), 20))

    transfer_all(left_pipette, temp_plate, thermocycler_plate, transfer_info_left3, mix_flag=1, mix_value=[5, 20])