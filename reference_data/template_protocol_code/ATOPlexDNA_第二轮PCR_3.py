from opentrons import protocol_api

metadata = {
    'protocolName': 'ATOPlexDNA_第二轮PCR_3',
    'author': 'MGIX',
    'description': 'Protocol for DNA fragmentation and recovery',
    'apiLevel': '2.14'
}


# Transfer function
def transfer_all(pipette, source_plate, dest_plate, transfer_info,mix_flag = 0,mix_value =[]):
    """
    General function for transferring samples.
    Args:
    - pipette: The pipette to use
    - source_plate: Source plate
    - dest_plate: Destination plate
    - transfer_info: List containing transfer information, each element is a tuple (source column name, destination column name, volume)
    - mix_value: Mixing information as a list [number of times, volume]
    """
    source_map = source_plate.columns_by_name()
    dest_map = dest_plate.columns_by_name()
    for source_col, dest_col, volume in transfer_info:
        temp_source = source_map.get(source_col)
        temp_dest = dest_map.get(dest_col)
        pipette.pick_up_tip()
        if mix_flag == 1:
            pipette.transfer(volume, temp_source, temp_dest, new_tip='never',mix_after=(mix_value[0], mix_value[1]))
        else:
            pipette.transfer(volume, temp_source, temp_dest, new_tip='never')
        pipette.drop_tip()


def calculate_liquid(sample_num,liquid_volume):
    liquid_list =[]
    if sample_num <8:
        per_sample_liquid = liquid_volume/sample_num
        temp_num = sample_num
        for ii in range(8):
            if temp_num > 0:
                liquid_list.append(per_sample_liquid)
                temp_num=temp_num-1
            else:
                liquid_list.append(0.0)
        return liquid_list
    else:
        count1 = int(sample_num/8)
        count2 = sample_num%8
        per_sample_liquid = liquid_volume/(8*(count1+1)+count2)

        for ii in range(8):
            if count2 > 0:
                liquid_list.append((count1+2)*per_sample_liquid)
                count2=count2-1
            else:
                liquid_list.append((count1+1)*per_sample_liquid)
        return liquid_list
def get_annealing_time(Amplicon_multiplicity):
    annealing_time = 0
    repet_times = 1
    if Amplicon_multiplicity<50:
        annealing_time = 30
        repet_times = 20
    elif Amplicon_multiplicity>=50 and Amplicon_multiplicity<200:
        annealing_time = 60
        repet_times = 19
    elif Amplicon_multiplicity>=200 and Amplicon_multiplicity<500:
        annealing_time = 60
        repet_times = 18
    elif Amplicon_multiplicity>=500 and Amplicon_multiplicity<1000:
        annealing_time = 120
        repet_times = 17
    elif Amplicon_multiplicity>=1000 and Amplicon_multiplicity<2000:
        annealing_time = 120
        repet_times = 16
    elif Amplicon_multiplicity>=2000 and Amplicon_multiplicity<5000:
        annealing_time = 180
        repet_times = 15
    elif Amplicon_multiplicity>=5000 and Amplicon_multiplicity<10000:
        annealing_time = 180
        repet_times = 14
    else:
        annealing_time = 240
        repet_times = 13
    return annealing_time,repet_times

def run(protocol: protocol_api.ProtocolContext):
    Sample_nums = input_sample_nums_1
    MIX_volume = input_PCR_MIX_2nd_volume_3
    
    Spike_in_Control_PCR_block = input_Spike_in_Control_PCR_block_1
    Amplicon_multiplicity = input_amplicon_multiplicity_1
    
    mix_per_sample = 14.5
    if Spike_in_Control_PCR_block:
        mix_per_sample = mix_per_sample +1
    annealing_time,repet_times = get_annealing_time(Amplicon_multiplicity)
    
    
    use_col_nums = int(Sample_nums/8)
    if Sample_nums%8 !=0:
        use_col_nums = use_col_nums+1
    # Load 20ul tip rack
    tips20_1 = protocol.load_labware('opentrons_96_tiprack_20ul', '2')
    tips20_2 = protocol.load_labware('opentrons_96_tiprack_20ul', '4')

    plate1 = protocol.load_labware('biorad_96_wellplate_200ul_pcr', '1')

    # Load temperature module and plate on temperature module
    temp_module = protocol.load_module('temperature module', '3')
    temp_plate = temp_module.load_labware('biorad_96_wellplate_200ul_pcr')

    # Load thermocycler module, ensure correct position number
    thermocycler_module = protocol.load_module('thermocycler', '7')
    thermocycler_plate = thermocycler_module.load_labware('biorad_96_wellplate_200ul_pcr')
    # Pipettes
    left_pipette = protocol.load_instrument('p20_multi_gen2', mount='left', tip_racks=[tips20_1,tips20_2])
    left_pipette.flow_rate.aspirate = 3.78  # Set aspirate speed to 50 μL/s
    left_pipette.flow_rate.dispense = 7.56 
    
    # Define liquids
    # labeling liquids in wells
    PCR_barcode = protocol.define_liquid(
        name="PCR_barcode",
        description="PCR_barcode",
        display_color="#00FF00",
    )
    PCR_mix = protocol.define_liquid(
        name="PCR_mix",
        description="PCR_mix",
        display_color="#00FF00",
    )
    PCR_Product = protocol.define_liquid(
        name="PCR_Product",
        description="PCR_Product",
        display_color="#00FF00",
    )
    # Load liquids
    all_plate_well = [f"{chr(65 + i)}{j}" for j in range(1, 12 + 1) for i in range(8)]
    liquid_list = calculate_liquid(Sample_nums,MIX_volume)
    for ii in range(8):
        well_name = all_plate_well[ii]
        temp_plate.wells_by_name()[well_name].load_liquid(liquid=PCR_mix, volume=liquid_list[ii])

    use_plate_well = all_plate_well[:Sample_nums]
    for well_name in use_plate_well:
        plate1.wells_by_name()[well_name].load_liquid(liquid=PCR_barcode, volume=30)
        thermocycler_plate.wells_by_name()[well_name].load_liquid(liquid=PCR_Product, volume=6.5)

    # Execute commands
    # 1. Set temperature module to 4°C and wait until temperature is reached
    temp_module.set_temperature(4)
    temp_module.await_temperature(4)

    thermocycler_module.open_lid()
    thermocycler_module.set_block_temperature(4)
    # Transfer
    transfer_times = use_col_nums
    transfer_info_left1 = []
    for ii in range(transfer_times):
        transfer_info_left1.append((str(ii + 1), str(ii + 1), mix_per_sample))
    transfer_all(left_pipette, plate1, thermocycler_plate, transfer_info_left1)

    transfer_info_left2 = []
    for ii in range(transfer_times):
        transfer_info_left2.append(('1',str(ii + 1), 4))

    transfer_all(left_pipette, temp_plate, thermocycler_plate, transfer_info_left2,mix_flag=1,mix_value=[5,10])

    thermocycler_module.close_lid()
    thermocycler_module.set_lid_temperature(105)
    steps1 = [{'temperature': 37, 'hold_time_minutes': 5},{'temperature': 95, 'hold_time_minutes': 10}]
    steps2 = [{'temperature': 95, 'hold_time_seconds': 20}, {'temperature': 62, 'hold_time_seconds': annealing_time},
              {'temperature': 58, 'hold_time_seconds': annealing_time},{'temperature': 72, 'hold_time_seconds': 30}]
    steps3 = [{'temperature': 4, 'hold_time_minutes': 1}]
    thermocycler_module.execute_profile(steps=steps1, repetitions=1,block_max_volume=25)
    thermocycler_module.execute_profile(steps=steps2, repetitions=repet_times,block_max_volume=25)
    thermocycler_module.execute_profile(steps=steps3, repetitions=1,block_max_volume=25)
    thermocycler_module.set_block_temperature(4, hold_time_seconds=None)
    thermocycler_module.deactivate_lid()
    thermocycler_module.open_lid()