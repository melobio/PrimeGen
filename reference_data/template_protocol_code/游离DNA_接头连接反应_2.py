from opentrons import protocol_api

metadata = {
    'protocolName': '游离DNA-2',
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
        per_sample_liquid = liquid_volume / (8 * (count1 + 0.5) + count2)

        for ii in range(8):
            if count2 > 0:
                liquid_list.append((count1 + 1.5) * per_sample_liquid)
                count2 = count2 - 1
            else:
                liquid_list.append((count1 + 1) * per_sample_liquid)
        return liquid_list

def run(protocol: protocol_api.ProtocolContext):
    Sample_nums = input_sample_nums_1
    MIX_volume = input_Adapter_ligation_volume_2
    use_col_nums = int(Sample_nums / 8)
    if Sample_nums % 8 != 0:
        use_col_nums = use_col_nums + 1
    # Load 20ul tip rack
    tips20 = protocol.load_labware('opentrons_96_tiprack_20ul', '6')

    # Load 200ul filter tip rack
    tips200_1 = protocol.load_labware('opentrons_96_tiprack_300ul', '2')
    tips200_2 = protocol.load_labware('opentrons_96_tiprack_300ul', '4')
    tips200_3 = protocol.load_labware('opentrons_96_tiprack_300ul', '5')

    plate1 = protocol.load_labware('biorad_96_wellplate_200ul_pcr', '1')

    # Load temperature module and load plate on temperature module
    temp_module = protocol.load_module('temperature module', '3')
    temp_plate = temp_module.load_labware('biorad_96_wellplate_200ul_pcr')

    # Load thermocycler module, ensure to specify a correct slot number
    thermocycler_module = protocol.load_module('thermocycler', '7')
    thermocycler_plate = thermocycler_module.load_labware('biorad_96_wellplate_200ul_pcr')

    # Pipettes
    left_pipette = protocol.load_instrument('p20_multi_gen2', mount='left', tip_racks=[tips20])
    right_pipette = protocol.load_instrument('p300_multi_gen2', mount='right', tip_racks=[tips200_1, tips200_2, tips200_3])
    left_pipette.flow_rate.aspirate = 3.78  # Set aspirate speed to 50 μL/s
    left_pipette.flow_rate.dispense = 7.56 
    right_pipette.flow_rate.aspirate = 46.43  # Set aspirate speed to 50 μL/s
    right_pipette.flow_rate.dispense = 92.86 
    # Define liquids
    Enzyme_digestion_solution_after = protocol.define_liquid(
        name="Enzyme digestion reaction solution",
        description="Enzyme digestion reaction solution",
        display_color="#00FF00",
    )
    Adapter_ligation_solution = protocol.define_liquid(
        name="Adapter ligation reaction solution",
        description="Adapter ligation reaction solution",
        display_color="#00FF00",
    )
    DNA_adapter = protocol.define_liquid(
        name="DNA_adapter",
        description="DNA_adapter",
        display_color="#00FF00",
    )
    # Load liquids
    all_plate_well = [f"{chr(65 + i)}{j}" for j in range(1, 12 + 1) for i in range(8)]
    if use_col_nums > 6:
        for well in temp_plate.columns_by_name()['1']:
            temp_plate.wells_by_name()[well.display_name].load_liquid(liquid=Adapter_ligation_solution, volume=162.5)
        liquid_list = calculate_liquid(Sample_nums - 48, MIX_volume - 162.5)
        for ii in range(8):
            temp_plate.wells_by_name()[all_plate_well[8:16][ii]].load_liquid(liquid=Adapter_ligation_solution, volume=liquid_list[ii])
    else:
        liquid_list = calculate_liquid(Sample_nums, MIX_volume)
        for ii in range(8):
            temp_plate.wells_by_name()[all_plate_well[ii]].load_liquid(liquid=Adapter_ligation_solution, volume=liquid_list[ii])
    
#     for well in temp_plate.columns_by_name()['1']:
#         temp_plate.wells_by_name()[well.display_name].load_liquid(liquid=Adapter_ligation_solution, volume=162.5)
#     for well in temp_plate.columns_by_name()['2']:
#         temp_plate.wells_by_name()[well.display_name].load_liquid(liquid=Adapter_ligation_solution, volume=162.5)

    use_plate_well = all_plate_well[:Sample_nums]
    for well_name in use_plate_well:
        thermocycler_plate.wells_by_name()[well_name].load_liquid(liquid=Enzyme_digestion_solution_after, volume=50)
        plate1.wells_by_name()[well_name].load_liquid(liquid=DNA_adapter, volume=30)
    # Execute instructions
    # 1. Set the temperature module to 4 degrees and pause to wait for the set temperature to be reached
    temp_module.set_temperature(4)
    temp_module.await_temperature(4)

    thermocycler_module.open_lid()
    thermocycler_module.set_block_temperature(4)
    # Transfer
    transfer_times = use_col_nums
    transfer_info_right1 = []
    for ii in range(transfer_times):
        if ii < 6:
            transfer_info_right1.append(('1', str(ii + 1), 25))
        else:
            transfer_info_right1.append(('2', str(ii + 1), 25))
    transfer_all(right_pipette, temp_plate, thermocycler_plate, transfer_info_right1)


    transfer_info_left1 = []
    for ii in range(transfer_times):
        transfer_info_left1.append((str(ii + 1), str(ii + 1), 5))

    transfer_all(left_pipette, plate1, thermocycler_plate, transfer_info_left1)

    # Mix
    for ii in range(transfer_times):
        temp_dest = thermocycler_plate.columns_by_name()[str(ii + 1)]
        right_pipette.pick_up_tip()
        right_pipette.mix(5, 50, temp_dest)
        right_pipette.drop_tip()


    thermocycler_module.close_lid()
    thermocycler_module.set_lid_temperature(105)
    thermocycler_module.execute_profile(steps=[
        {'temperature': 23, 'hold_time_minutes': 20},
        {'temperature': 4, 'hold_time_minutes': 1}
    ], repetitions=1, block_max_volume=80)
    thermocycler_module.set_block_temperature(4, hold_time_seconds=None)
    thermocycler_module.deactivate_lid()
    thermocycler_module.open_lid()