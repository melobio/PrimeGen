from opentrons import protocol_api

metadata = {
    'protocolName': 'MGIEasy fast-5',
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
        
def get_PCR_repet(temp_DNA,DNA_ng):
    temp_time = 0 
    DNA_list = [1,5,10,25,50,100,200,500,1000]
    ng300_list = [11,9,8,7,6,5,3,3,3]
    ng600_list = [13,10,9,8,7,6,5,4,3]
    ng1000_list = [11,11,10,9,8,7,6,5,4]
    for ii in range(len(DNA_list)):
        if temp_DNA ==DNA_list[ii] and DNA_ng==300:
            temp_time = ng300_list[ii]
        elif temp_DNA ==DNA_list[ii] and DNA_ng==600:
            temp_time = ng600_list[ii]
        elif temp_DNA ==DNA_list[ii] and DNA_ng==1000:
            temp_time = ng1000_list[ii]
           
    return temp_time

def run(protocol: protocol_api.ProtocolContext):
    Sample_nums = input_sample_nums_1
    PCR_MIX_volume = input_PCR_MIX_volume_5
    UDB_MIX_volume = input_UDB_MIX_volume_5
    DNA_quality = input_DNA_quality_5
    DNA_NG = input_DNA_NG_5
    repet_times = get_PCR_repet(DNA_quality,DNA_NG)
    use_col_nums = int(Sample_nums/8)
    if Sample_nums%8 !=0:
        use_col_nums = use_col_nums+1
    # 加载20ul吸头架
    tips20 = protocol.load_labware('opentrons_96_tiprack_20ul', '1')
    # 加载200ul滤芯吸头架
    tips200_1 = protocol.load_labware('opentrons_96_tiprack_300ul', '2')
    # 加载温度模块和在温度模块上加载
    temp_module = protocol.load_module('temperature module', '3')
    temp_plate = temp_module.load_labware('biorad_96_wellplate_200ul_pcr')
    # 加载热循环仪模块，确保指定一个正确的位置编号
    thermocycler_module = protocol.load_module('thermocycler', '7')
    thermocycler_plate = thermocycler_module.load_labware('biorad_96_wellplate_200ul_pcr')

    # Pipettes
    left_pipette = protocol.load_instrument('p20_multi_gen2', mount='left', tip_racks=[tips20])
    right_pipette = protocol.load_instrument('p300_multi_gen2', mount='right', tip_racks=[tips200_1])
    left_pipette.flow_rate.aspirate = 3.78  # Set aspirate speed to 50 μL/s
    left_pipette.flow_rate.dispense = 7.56 
    right_pipette.flow_rate.aspirate = 46.43  # Set aspirate speed to 50 μL/s
    right_pipette.flow_rate.dispense = 92.86 
    
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
    UDB_MIX = protocol.define_liquid(
        name="UDB_MIX",
        description="UDB_MIX",
        display_color="#00FF00",
    )
    # 加载液体
    all_plate_well = [f"{chr(65 + i)}{j}" for j in range(1, 12 + 1) for i in range(8)]
    if use_col_nums>6:
        for well in temp_plate.columns_by_name()['1']:
            temp_plate.wells_by_name()[well.display_name].load_liquid(liquid=PCR_MIX, volume=162.5)
        liquid_list = calculate_liquid(Sample_nums-48,PCR_MIX_volume-162.5)
        for ii in range(8):
            temp_plate.wells_by_name()[all_plate_well[8:16][ii]].load_liquid(liquid=PCR_MIX, volume=liquid_list[ii])
    else:
        liquid_list = calculate_liquid(Sample_nums,PCR_MIX_volume)
        for ii in range(8):
            temp_plate.wells_by_name()[all_plate_well[ii]].load_liquid(liquid=PCR_MIX, volume=liquid_list[ii])
    
#     for well in temp_plate.columns_by_name()['1']:
#         temp_plate.wells_by_name()[well.display_name].load_liquid(liquid=PCR_MIX, volume=162.5)
#     for well in temp_plate.columns_by_name()['2']:
#         temp_plate.wells_by_name()[well.display_name].load_liquid(liquid=PCR_MIX, volume=162.5)
    liquid_list = calculate_liquid(Sample_nums,UDB_MIX_volume)
    for ii in range(8):
        well_name = all_plate_well[16:24][ii]
        temp_plate.wells_by_name()[well_name].load_liquid(liquid=UDB_MIX, volume=liquid_list[ii]) 
        
#     for well in temp_plate.columns_by_name()['3']:
#         temp_plate.wells_by_name()[well.display_name].load_liquid(liquid=UDB_MIX, volume=78)
        
    use_plate_well = all_plate_well[:Sample_nums]
    for well_name in use_plate_well: 
        thermocycler_plate.wells_by_name()[well_name].load_liquid(liquid=sample, volume=19)

    #执行指令
    thermocycler_module.deactivate_block()
    thermocycler_module.deactivate_lid()
    thermocycler_module.open_lid()
    # 移液
    transfer_times =use_col_nums
    transfer_info_left1 = []
    for ii in range(transfer_times):
        transfer_info_left1.append(('3', str(ii + 1), 6))

    transfer_all(left_pipette, temp_plate, thermocycler_plate, transfer_info_left1)

    transfer_info_right1 = []
    for ii in range(transfer_times):
        if ii<6:
            transfer_info_right1.append(('1',str(ii+1),25))
        else:
            transfer_info_right1.append(('2', str(ii + 1), 25))
    
    transfer_all(right_pipette, temp_plate, thermocycler_plate, transfer_info_right1,mix_flag=1,mix_value=[5,25])

    thermocycler_module.close_lid()
    thermocycler_module.set_lid_temperature(105)
    steps1 = [{'temperature': 95, 'hold_time_minutes': 3}]
    steps2 = [{'temperature': 98, 'hold_time_seconds': 20},{'temperature': 60, 'hold_time_seconds': 15},{'temperature': 72, 'hold_time_seconds': 30}]
    steps3 = [{'temperature': 72, 'hold_time_minutes': 10},{'temperature': 4, 'hold_time_minutes': 1}]
    thermocycler_module.execute_profile(steps=steps1, repetitions=1,block_max_volume=50)
    thermocycler_module.execute_profile(steps=steps2, repetitions=repet_times,block_max_volume=50)
    thermocycler_module.execute_profile(steps=steps3, repetitions=1,block_max_volume=50)
    thermocycler_module.set_block_temperature(4, hold_time_seconds=None)
    thermocycler_module.deactivate_lid()
    thermocycler_module.open_lid()

