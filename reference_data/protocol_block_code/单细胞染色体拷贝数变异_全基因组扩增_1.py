from opentrons import protocol_api

metadata = {
    'protocolName': '单细胞染色体拷贝数变异_全基因组扩增_1',
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
    - pipette: 使用的移液器
    - source_plate: 源板
    - dest_plate: 目标板
    - transfer_info: 包含移液信息的列表，每个元素为一个元组 (源列名, 目标列名, 体积)
    - mix_value: 混液信息，是一个列表，[次数，体积]
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
    Sample_nums = 96
    Cell_lysis_master_mix_volume = 2808
    Pre_PCR_master_mix_volume = 2808
    Post_PCR_master_mix_volume = 2808
    Nuclease_free_water_volume = 2808

    use_col_nums = int(Sample_nums / 8)
    if Sample_nums % 8 != 0:
        use_col_nums = use_col_nums + 1
    # Load 200ul tip rack
    tips20_1 = protocol.load_labware('opentrons_96_tiprack_20ul', '1')
    tips20_2 = protocol.load_labware('opentrons_96_tiprack_20ul', '2')

    tips200_1 = protocol.load_labware('opentrons_96_tiprack_300ul', '4')
    tips200_2 = protocol.load_labware('opentrons_96_tiprack_300ul', '5')
    temp_module = protocol.load_module('temperature module', '3')
    temp_plate = temp_module.load_labware('biorad_96_wellplate_200ul_pcr')
    # Load thermocycler module, ensure to specify the correct position number
    thermocycler_module = protocol.load_module('thermocycler', '7')
    thermocycler_plate = thermocycler_module.load_labware('biorad_96_wellplate_200ul_pcr')
    # Pipettes
    left_pipette = protocol.load_instrument('p20_multi_gen2', mount='left', tip_racks=[tips20_1,tips20_2])
    right_pipette = protocol.load_instrument('p300_multi_gen2', mount='right', tip_racks=[tips200_1,tips200_2])
    right_pipette.flow_rate.aspirate = 46.43  # Set aspirate speed to 50 μL/s
    right_pipette.flow_rate.dispense = 92.86
    left_pipette.flow_rate.aspirate = 3.78  # Set aspirate speed to 50 μL/s
    left_pipette.flow_rate.dispense = 7.56
    # Define liquids
    # labeling liquids in wells
    sample = protocol.define_liquid(
        name="sample",
        description="sample",
        display_color="#00FF00",
    )
    Cell_lysis_master_mix = protocol.define_liquid(
        name="Cell_lysis_master_mix",
        description="Cell_lysis_master_mix",
        display_color="#00FF00",
    )
    Pre_PCR_master_mix = protocol.define_liquid(
        name="Pre_PCR_master_mix",
        description="Pre_PCR_master_mix",
        display_color="#00FF00",
    )
    Post_PCR_master_mix = protocol.define_liquid(
        name="Post_PCR_master_mix",
        description="Post_PCR_master_mix",
        display_color="#00FF00",
    )
    Nuclease_free_water= protocol.define_liquid(
        name="Nuclease_free_water",
        description="Nuclease_free_water",
        display_color="#00FF00",
    )
    # Load liquids
    all_plate_well = [f"{chr(65 + i)}{j}" for j in range(1, 12 + 1) for i in range(8)]
    liquid_list1 = calculate_liquid(Sample_nums, Cell_lysis_master_mix_volume)
    liquid_list2 = calculate_liquid(Sample_nums, Pre_PCR_master_mix_volume)
    for ii in range(8):
        temp_plate.wells_by_name()[all_plate_well[ii]].load_liquid(liquid=Cell_lysis_master_mix, volume=liquid_list1[ii])
        temp_plate.wells_by_name()[all_plate_well[8:16][ii]].load_liquid(liquid=Pre_PCR_master_mix,
                                                                   volume=liquid_list2[ii])

    if use_col_nums > 6:
        for ii in range(8):
            temp_plate.wells_by_name()[all_plate_well[16:24][ii]].load_liquid(liquid=Post_PCR_master_mix, volume=167.7)
        liquid_list3 = calculate_liquid(Sample_nums - 48, Post_PCR_master_mix_volume - 167.7)
        for ii in range(8):
            temp_plate.wells_by_name()[all_plate_well[24:32][ii]].load_liquid(liquid=Post_PCR_master_mix, volume=liquid_list3[ii])
    else:
        liquid_list3 = calculate_liquid(Sample_nums, Post_PCR_master_mix_volume)
        for ii in range(8):
            temp_plate.wells_by_name()[all_plate_well[16:24][ii]].load_liquid(liquid=Post_PCR_master_mix, volume=liquid_list3[ii])

    for well in temp_plate.columns_by_name()['5']:
        temp_plate.wells_by_name()[well.display_name].load_liquid(liquid=Nuclease_free_water, volume=148.2)
    for well in temp_plate.columns_by_name()['6']:
        temp_plate.wells_by_name()[well.display_name].load_liquid(liquid=Nuclease_free_water, volume=148.2)
    for well in temp_plate.columns_by_name()['7']:
        temp_plate.wells_by_name()[well.display_name].load_liquid(liquid=Nuclease_free_water, volume=148.2)

    use_plate_well = all_plate_well[:Sample_nums]
    for well_name in use_plate_well:
        thermocycler_plate.wells_by_name()[well_name].load_liquid(liquid=sample, volume=4)

    # Execute instructions
    # 1. Set the temperature module to 4°C and pause until the set temperature is reached
    temp_module.set_temperature(4)
    temp_module.await_temperature(4)
    # Open thermocycler
    thermocycler_module.open_lid()
    thermocycler_module.set_block_temperature(4)
    # Liquid transfer
    transfer_times = use_col_nums
    transfer_info_left1 = []
    for ii in range(transfer_times):
        transfer_info_left1.append(('1', str(ii + 1), 5))
    transfer_all(left_pipette, temp_plate, thermocycler_plate, transfer_info_left1, mix_flag=1, mix_value=[10, 5])

    thermocycler_module.close_lid()
    thermocycler_module.set_lid_temperature(105)
    steps1 = [{'temperature': 75, 'hold_time_minutes': 10},{'temperature': 95, 'hold_time_minutes': 4}, {'temperature': 4, 'hold_time_minutes': 1}]
    thermocycler_module.execute_profile(steps=steps1, repetitions=1, block_max_volume=9)
    thermocycler_module.set_block_temperature(4, hold_time_seconds=None)
    thermocycler_module.deactivate_lid()
    thermocycler_module.open_lid()

    transfer_info_left2 = []
    for ii in range(transfer_times):
        transfer_info_left2.append(('2', str(ii + 1), 5))
    transfer_all(left_pipette, temp_plate, thermocycler_plate, transfer_info_left2, mix_flag=1, mix_value=[5, 10])

    thermocycler_module.close_lid()
    thermocycler_module.set_lid_temperature(105)
    steps1 = [{'temperature': 95, 'hold_time_minutes': 2}]
    steps2 = [{'temperature': 95, 'hold_time_seconds': 15}, {'temperature': 15, 'hold_time_seconds': 50},
              {'temperature': 25, 'hold_time_seconds': 40},{'temperature': 35, 'hold_time_seconds': 30},
              {'temperature': 65, 'hold_time_seconds': 40},{'temperature': 75, 'hold_time_seconds': 40}]
    steps3 = [{'temperature': 4, 'hold_time_minutes': 1}]
    thermocycler_module.execute_profile(steps=steps1, repetitions=1, block_max_volume=14)
    thermocycler_module.execute_profile(steps=steps2, repetitions=12, block_max_volume=14)
    thermocycler_module.execute_profile(steps=steps3, repetitions=1, block_max_volume=14)
    thermocycler_module.set_block_temperature(4, hold_time_seconds=None)
    thermocycler_module.deactivate_lid()
    thermocycler_module.open_lid()

    transfer_info_right1 = []
    for ii in range(transfer_times):
        if ii < 6:
            transfer_info_right1.append(('3', str(ii + 1), 25.8))
        else:
            transfer_info_right1.append(('4', str(ii + 1), 25.8))
    transfer_all(right_pipette, temp_plate, thermocycler_plate, transfer_info_right1, mix_flag=1, mix_value=[5, 20])

    transfer_info_right2 = []
    for ii in range(transfer_times):
        if ii < 4:
            transfer_info_right2.append(('5', str(ii + 1), 34.2))
        elif ii>=4 and ii<8:
            transfer_info_right2.append(('6', str(ii + 1), 34.2))
        else:
            transfer_info_right2.append(('7', str(ii + 1), 34.2))
    transfer_all(right_pipette, temp_plate, thermocycler_plate, transfer_info_right2, mix_flag=1, mix_value=[5, 50])


    thermocycler_module.close_lid()
    thermocycler_module.set_lid_temperature(105)
    steps1 = [{'temperature': 95, 'hold_time_minutes': 2}]
    steps2 = [{'temperature': 95, 'hold_time_seconds': 15}, {'temperature': 65, 'hold_time_seconds': 60},
              {'temperature': 75, 'hold_time_seconds': 60}]
    steps3 = [{'temperature': 4, 'hold_time_minutes': 1}]
    thermocycler_module.execute_profile(steps=steps1, repetitions=1, block_max_volume=74)
    thermocycler_module.execute_profile(steps=steps2, repetitions=14, block_max_volume=74)
    thermocycler_module.execute_profile(steps=steps3, repetitions=1, block_max_volume=74)
    thermocycler_module.set_block_temperature(4, hold_time_seconds=None)
    thermocycler_module.deactivate_lid()
    thermocycler_module.open_lid()
