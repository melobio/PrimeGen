from opentrons import protocol_api

metadata = {
    'protocolName': 'MGIEasy_fastRNA_RNA二链合成_4',
    'author': 'MGIX',
    'description': 'Protocol for DNA fragmentation and recovery',
    'apiLevel': '2.14'
}


# Liquid transfer function
def transfer_all(pipette, source_plate, dest_plate, transfer_info, mix_flag=0, mix_value=[]):
    """
    General function for transferring samples.
    Args:
    - pipette: The pipette used
    - source_plate: Source plate
    - dest_plate: Destination plate
    - transfer_info: A list containing transfer information, where each element is a tuple (source column name, destination column name, volume)
    - mix_value: Mixing information, which is a list [number of times, volume]
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


def calculate_liquid(sample_num, liquid_volume,cal_type=1):
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
        if cal_type ==1:
            count1 = int(sample_num / 8)
            count2 = sample_num % 8
            per_sample_liquid = liquid_volume / (8 * (count1 + 1) + count2)

            for ii in range(8):
                if count2 > 0:
                    liquid_list.append((count1 + 2) * per_sample_liquid)
                    count2 = count2 - 1
                else:
                    liquid_list.append((count1 + 1) * per_sample_liquid)
        else:
            count1 = int(sample_num / 8)
            count2 = sample_num % 8
            per_sample_liquid = liquid_volume / (8 * (count1 + 0.5) + count2)

            for ii in range(8):
                if count2 > 0:
                    liquid_list.append((count1 + 1.5) * per_sample_liquid)
                    count2 = count2 - 1
                else:
                    liquid_list.append((count1 + 0.5) * per_sample_liquid)
        return liquid_list

def run(protocol: protocol_api.ProtocolContext):
    Sample_nums = user_sample_nums
    Second_strand_master_mix_volume = input_Second_strand_master_mix_volume
    use_col_nums = int(Sample_nums / 8)
    if Sample_nums % 8 != 0:
        use_col_nums = use_col_nums + 1
    # Load 20ul tip rack
    tips200_1 = protocol.load_labware('opentrons_96_tiprack_300ul', '1')
    temp_module = protocol.load_module('temperature module', '3')
    temp_plate = temp_module.load_labware('biorad_96_wellplate_200ul_pcr')
    # Load thermocycler module, ensure correct position number
    thermocycler_module = protocol.load_module('thermocycler', '7')
    thermocycler_plate = thermocycler_module.load_labware('biorad_96_wellplate_200ul_pcr')
    # Pipettes
    right_pipette = protocol.load_instrument('p300_multi_gen2', mount='right', tip_racks=[tips200_1])
    right_pipette.flow_rate.aspirate = 46.43  # Set aspirate speed to 50 μL/s
    right_pipette.flow_rate.dispense = 92.86

    # Define liquids
    # labeling liquids in wells
    RNA_sample = protocol.define_liquid(
        name="RNA_sample",
        description="RNA_sample",
        display_color="#00FF00",
    )
    Second_strand_master_mix = protocol.define_liquid(
        name="Second_strand_master_mix",
        description="Second_strand_master_mix",
        display_color="#00FF00",
    )
    # Load liquids
    all_plate_well = [f"{chr(65 + i)}{j}" for j in range(1, 12 + 1) for i in range(8)]
    if use_col_nums > 6:
        half_volume = (Second_strand_master_mix_volume / (Sample_nums + 8)) * 6.5
        for ii in range(8):
            temp_plate.wells_by_name()[all_plate_well[ii]].load_liquid(liquid=Second_strand_master_mix, volume=half_volume)
        liquid_list = calculate_liquid(Sample_nums - 48, Second_strand_master_mix_volume - half_volume * 8, cal_type=2)
        for ii in range(8):
            temp_plate.wells_by_name()[all_plate_well[8:16][ii]].load_liquid(liquid=Second_strand_master_mix, volume=liquid_list[ii])
    else:
        liquid_list = calculate_liquid(Sample_nums, Second_strand_master_mix_volume)
        for ii in range(8):
            temp_plate.wells_by_name()[all_plate_well[ii]].load_liquid(liquid=Second_strand_master_mix, volume=liquid_list[ii])

    use_plate_well = all_plate_well[:Sample_nums]
    for well_name in use_plate_well:
        thermocycler_plate.wells_by_name()[well_name].load_liquid(liquid=RNA_sample, volume=20)
        # Execute commands
    # 1. Set temperature module to 4°C and wait until target temperature is reached
    temp_module.set_temperature(4)
    temp_module.await_temperature(4)
    # Start thermocycler
    thermocycler_module.open_lid()
    thermocycler_module.set_block_temperature(4)
    # Transfer
    transfer_times = use_col_nums
    transfer_info_right1 = []
    for ii in range(transfer_times):
        if ii < 6:
            transfer_info_right1.append(('1', str(ii + 1), 30))
        else:
            transfer_info_right1.append(('2', str(ii + 1), 30))
    transfer_all(right_pipette, temp_plate, thermocycler_plate, transfer_info_right1, mix_flag=1, mix_value=[10, 15])

    thermocycler_module.close_lid()
    thermocycler_module.set_lid_temperature(70)
    steps1 = [{'temperature': 16, 'hold_time_minutes': 30},{'temperature': 65, 'hold_time_minutes': 15},
              {'temperature': 4, 'hold_time_minutes': 2}]
    thermocycler_module.execute_profile(steps=steps1, repetitions=1, block_max_volume=50)
    thermocycler_module.set_block_temperature(4, hold_time_seconds=None)
    thermocycler_module.deactivate_lid()
    thermocycler_module.open_lid()