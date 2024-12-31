from opentrons import protocol_api
##
## This code block is suitable for PCR amplification reactions in protocols such as ATOPlexHIV_PCR, Multiplex Amplification_1.py, ATOPlexMPVX_FirstRoundPCRPOOL_2.py, ATOPlexMPVX_FirstRoundPCRPOOL_1.py, etc.
##

metadata = {
    'protocolName': 'ATOPlexMPVX-1',
    'author': 'MGIX',
    'description': 'Protocol for DNA fragmentation and recovery',
    'apiLevel': '2.14'
}

# Liquid handling function
def transfer_all(pipette, source_plate, dest_plate, transfer_info, mix_flag=0, mix_value=[]):
    """
    General function for transferring samples.
    Args:
    - pipette: Pipette used
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

        
def calculate_liquid(sample_num, liquid_volume, cal_type=1):
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
        if cal_type == 1:
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
    MIX_volume = input_PCR_MIX_volume
    Single_sample_volume = input_Single_sample_volume
    thermocycler_param = input_thermocycler_param
    
    transfer_volume_1 = 10
    transfer_volume_2 = MIX_volume / (Sample_nums + 8)
    

#     [{'temper_time':[(37,300),(95,300)],'repetition':1},
#                     {'temper_time':[(95,20),(62,60),(58,60),(72,20)],'repetition':37},
#                     {'temper_time':[(4,60)],'repetition':1}]
    use_col_nums = int(Sample_nums / 8)
    if Sample_nums % 8 != 0:
        use_col_nums = use_col_nums + 1

    temp_module = protocol.load_module('temperature module', '3')
    temp_plate = temp_module.load_labware('biorad_96_wellplate_200ul_pcr')
    plate1 = protocol.load_labware('biorad_96_wellplate_200ul_pcr', '1')

    # Load thermocycler module, ensure a correct slot number is specified
    thermocycler_module = protocol.load_module('thermocycler', '7')
    thermocycler_plate = thermocycler_module.load_labware('biorad_96_wellplate_200ul_pcr')
    
    if transfer_volume_2 > 20:
        # Load 20μL tip rack
        tips20_1 = protocol.load_labware('opentrons_96_tiprack_20ul', '2')
        tips200_1 = protocol.load_labware('opentrons_96_tiprack_300ul', '4')
        left_pipette = protocol.load_instrument('p20_multi_gen2', mount='left', tip_racks=[tips20_1])
        right_pipette = protocol.load_instrument('p300_multi_gen2', mount='right', tip_racks=[tips200_1])
        left_pipette.flow_rate.aspirate = 3.78  # Set aspirate speed to 50 μL/s
        left_pipette.flow_rate.dispense = 7.56 
        right_pipette.flow_rate.aspirate = 46.43  # Set aspirate speed to 50 μL/s
        right_pipette.flow_rate.dispense = 92.86 
    else:
        # Load 20μL tip rack
        tips20_1 = protocol.load_labware('opentrons_96_tiprack_20ul', '2')
        tips20_2 = protocol.load_labware('opentrons_96_tiprack_20ul', '4')
        # Pipettes
        left_pipette = protocol.load_instrument('p20_multi_gen2', mount='left', tip_racks=[tips20_1, tips20_2])
        left_pipette.flow_rate.aspirate = 3.78  # Set aspirate speed to 50 μL/s
        left_pipette.flow_rate.dispense = 7.56 
    # Define liquids
    # Labeling liquids in wells
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
    all_plate_well = [f"{chr(65 + i)}{j}" for j in range(1, 12 + 1) for i in range(8)]
    if transfer_volume_2 > 20:
        
        if use_col_nums > 6:
            half_volume = (MIX_volume / (Sample_nums + 8)) * 6.5
            for ii in range(8):
                temp_plate.wells_by_name()[all_plate_well[ii]].load_liquid(liquid=PCR_MIX, volume=half_volume)
            liquid_list = calculate_liquid(Sample_nums - 48, MIX_volume - half_volume * 8, cal_type=2)
            for ii in range(8):
                temp_plate.wells_by_name()[all_plate_well[8:16][ii]].load_liquid(liquid=PCR_MIX, volume=liquid_list[ii])
        else:
            liquid_list = calculate_liquid(Sample_nums, MIX_volume)
            for ii in range(8):
                temp_plate.wells_by_name()[all_plate_well[ii]].load_liquid(liquid=PCR_MIX, volume=liquid_list[ii])
        
    else:
        liquid_list = calculate_liquid(Sample_nums, MIX_volume)
        for ii in range(8):
            well_name = all_plate_well[ii]
            temp_plate.wells_by_name()[well_name].load_liquid(liquid=PCR_MIX, volume=liquid_list[ii])

    use_plate_well = all_plate_well[:Sample_nums]
    for well_name in use_plate_well:
        plate1.wells_by_name()[well_name].load_liquid(liquid=sample, volume=Single_sample_volume)

    # Execute commands
    # 1. Set the temperature module to 4°C and pause until the target temperature is reached
    temp_module.set_temperature(4)
    temp_module.await_temperature(4)
    # Start thermocycler
    thermocycler_module.open_lid()
    thermocycler_module.set_block_temperature(4)
    # Liquid handling
    transfer_times = use_col_nums
    transfer_info_left1 = []
    for ii in range(transfer_times):
        transfer_info_left1.append((str(ii + 1), str(ii + 1), transfer_volume_1))

    transfer_all(left_pipette, plate1, thermocycler_plate, transfer_info_left1)

    transfer_info_2 = []
    if transfer_volume_2 > 20:
        for ii in range(transfer_times):
            if ii < 6:
                transfer_info_2.append(('1', str(ii + 1), transfer_volume_2))
            else:
                transfer_info_2.append(('2', str(ii + 1), transfer_volume_2))
        transfer_all(right_pipette, temp_plate, thermocycler_plate, transfer_info_2, mix_flag=1, mix_value=[10, 20])

    else:
        for ii in range(transfer_times):
            transfer_info_2.append(('1', str(ii + 1), transfer_volume_1))
        transfer_all(left_pipette, temp_plate, thermocycler_plate, transfer_info_2, mix_flag=1, mix_value=[5, 10])

    thermocycler_module.close_lid()
    thermocycler_module.set_lid_temperature(105)
    block_max_volume = 25
    if transfer_volume_2 > 20:
        block_max_volume = 50
    
    for param in thermocycler_param:
        temper_time = param['param']
        repetitions = param['repetitions']
        temp_step = []
        for tem in temper_time:
            temp_step.append({'temperature': tem[0], 'hold_time_seconds': tem[1]})
            
        thermocycler_module.execute_profile(steps=temp_step, repetitions=repetitions, block_max_volume=block_max_volume)
    thermocycler_module.set_block_temperature(4, hold_time_seconds=None)
    thermocycler_module.deactivate_lid()
    thermocycler_module.open_lid()
