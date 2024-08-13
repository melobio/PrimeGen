from opentrons import protocol_api

metadata = {
    'protocolName': 'ATOPlexRNA-2',
    'author': 'MGIX',
    'description': 'Protocol for DNA fragmentation and recovery',
    'apiLevel': '2.14'
}


#移液函数
def transfer_all(pipette, source_plate, dest_plate, transfer_info,mix_flag = 0,mix_value =[]):
    """
    传输样本的通用函数。
    Args:
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
            pipette.transfer(volume, temp_source, temp_dest, new_tip='never',mix_after=(mix_value[0], mix_value[1]))
        else:
            pipette.transfer(volume, temp_source, temp_dest, new_tip='never')
        pipette.drop_tip()
        
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

def run(protocol: protocol_api.ProtocolContext):
    Sample_nums = 96
    MIX_volume = 1560
    use_col_nums = int(Sample_nums/8)
    if Sample_nums%8 !=0:
        use_col_nums = use_col_nums+1
        
    Spike_in_Control_PCR_block = False
    Amplicon_multiplicity = 50
    mix_per_sample = 15
    
    if Spike_in_Control_PCR_block:
        mix_per_sample = mix_per_sample +1
    annealing_time,repet_times = get_annealing_time(Amplicon_multiplicity)
        
    # 加载20ul吸头架
    tips20_1 = protocol.load_labware('opentrons_96_tiprack_20ul', '2')
    tips20_2 = protocol.load_labware('opentrons_96_tiprack_20ul', '4')

    temp_module = protocol.load_module('temperature module', '3')
    temp_plate = temp_module.load_labware('biorad_96_wellplate_200ul_pcr')
    plate1 = protocol.load_labware('biorad_96_wellplate_200ul_pcr', '1')
    # 加载热循环仪模块，确保指定一个正确的位置编号
    thermocycler_module = protocol.load_module('thermocycler', '7')
    thermocycler_plate = thermocycler_module.load_labware('biorad_96_wellplate_200ul_pcr')
    # Pipettes
    left_pipette = protocol.load_instrument('p20_multi_gen2', mount='left', tip_racks=[tips20_1,tips20_2])
    left_pipette.flow_rate.aspirate = 3.78  # Set aspirate speed to 50 μL/s
    left_pipette.flow_rate.dispense = 7.56 
    # 定义液体
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
    # 加载液体
    all_plate_well = [f"{chr(65 + i)}{j}" for j in range(1, 12 + 1) for i in range(8)]
    
    liquid_list = calculate_liquid(Sample_nums,MIX_volume)
    for ii in range(8):
        well_name = all_plate_well[ii]
        temp_plate.wells_by_name()[well_name].load_liquid(liquid=PCR_MIX, volume=liquid_list[ii])
        
#     for well in temp_plate.columns_by_name()['1']:
#         temp_plate.wells_by_name()[well.display_name].load_liquid(liquid=PCR_MIX, volume=195)

    use_plate_well = all_plate_well[:Sample_nums]
    for well_name in use_plate_well: 
        plate1.wells_by_name()[well_name].load_liquid(liquid=sample, volume=20)

    #执行指令
    # 1. 将温控模块设置为4度并暂停等待达到设定温度
    temp_module.set_temperature(4)
    temp_module.await_temperature(4)
    # 开启热循环
    thermocycler_module.open_lid()
    thermocycler_module.set_block_temperature(4)
    # 移液
    transfer_times = use_col_nums
    transfer_info_left1 = []
    for ii in range(transfer_times):
        transfer_info_left1.append((str(ii + 1), str(ii + 1), 10))

    transfer_all(left_pipette, plate1, thermocycler_plate, transfer_info_left1)

    transfer_info_left2 = []
    for ii in range(transfer_times):
        transfer_info_left2.append(('1', str(ii + 1), mix_per_sample))

    transfer_all(left_pipette, temp_plate, thermocycler_plate, transfer_info_left2, mix_flag=1, mix_value=[5, 10])

    thermocycler_module.close_lid()
    thermocycler_module.set_lid_temperature(105)
    steps1 = [{'temperature': 37, 'hold_time_minutes': 5},{'temperature': 95, 'hold_time_minutes': 10}]
    steps2 = [{'temperature': 95, 'hold_time_seconds': 20}, {'temperature': 62, 'hold_time_seconds': annealing_time},
              {'temperature': 58, 'hold_time_seconds': annealing_time},{'temperature': 72, 'hold_time_seconds': 30}]
    steps3 = [{'temperature': 4, 'hold_time_minutes': 1}]
    thermocycler_module.execute_profile(steps=steps1, repetitions=1,block_max_volume=25)
    thermocycler_module.execute_profile(steps=steps2, repetitions=10,block_max_volume=25)
    thermocycler_module.execute_profile(steps=steps3, repetitions=1,block_max_volume=25)
    thermocycler_module.set_block_temperature(4, hold_time_seconds=None)
    thermocycler_module.deactivate_lid()
    thermocycler_module.open_lid()
