from opentrons import protocol_api

metadata = {
    'protocolName': 'MGIEasy_fast_酶切打断_1',
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
def get_Enzyme_digestion_time(temp_DNA):
    temp_time = 0 
    temp_type = ''
    DNA_list = [1,5,10,25,50,100,200,500,1000]
    time_list = [22,20,18,16,16,13,12,12,12]
    type_list = ['磁珠单选','磁珠单选','磁珠单选','磁珠单选','磁珠单选','磁珠单选','磁珠单选','磁珠双选','磁珠双选']
    for ii in range(len(DNA_list)):
        if temp_DNA ==DNA_list[ii]:
            temp_time = time_list[ii]
            temp_type = type_list[ii]
            
    return temp_time,temp_type

def run(protocol: protocol_api.ProtocolContext):
    Sample_nums = 96
    MIX_volume = 1560
    DNA_quality = 50
    Enzyme_digestion_time,temp_type = get_Enzyme_digestion_time(DNA_quality)
    use_col_nums = int(Sample_nums/8)
    if Sample_nums%8 !=0:
        use_col_nums = use_col_nums+1
    # 加载20ul吸头架
    tips20 = protocol.load_labware('opentrons_96_tiprack_20ul', '4')
    # 加载200ul滤芯吸头架
    tips200_1 = protocol.load_labware('opentrons_96_tiprack_300ul', '2')
    plate1 = protocol.load_labware('biorad_96_wellplate_200ul_pcr', '1')

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
    # 定义液体
    # labeling liquids in wells
    sample = protocol.define_liquid(
        name="sample",
        description="sample",
        display_color="#00FF00",
    )
    Fast_FS_Enzyme = protocol.define_liquid(
        name="Fast_FS_Enzyme",
        description="Fast_FS_Enzyme",
        display_color="#00FF00",
    )
    # 加载液体
    all_plate_well = [f"{chr(65 + i)}{j}" for j in range(1, 12 + 1) for i in range(8)]
    liquid_list = calculate_liquid(Sample_nums,MIX_volume)
    for ii in range(8):
        well_name = all_plate_well[ii]
        temp_plate.wells_by_name()[well_name].load_liquid(liquid=Fast_FS_Enzyme, volume=liquid_list[ii])
        
#     for well in temp_plate.columns_by_name()['1']:
#         temp_plate.wells_by_name()[well.display_name].load_liquid(liquid=Fast_FS_Enzyme, volume=195)

    use_plate_well = all_plate_well[:Sample_nums]
    for well_name in use_plate_well: 
        plate1.wells_by_name()[well_name].load_liquid(liquid=sample, volume=50) 
    #执行指令
    #开启热循环
    thermocycler_module.open_lid()
    thermocycler_module.set_block_temperature(4)
    #移液
    transfer_times = use_col_nums
    transfer_info_right1 = []
    for ii in range(transfer_times):
        transfer_info_right1.append((str(ii+1),str(ii+1),45))

    transfer_all(right_pipette, plate1, thermocycler_plate, transfer_info_right1)


    transfer_info_left1 = []
    for ii in range(transfer_times):
        transfer_info_left1.append(('1', str(ii + 1), 15))

    transfer_all(left_pipette, temp_plate, thermocycler_plate, transfer_info_left1,mix_flag=1,mix_value=[5,10])

    thermocycler_module.close_lid()
    thermocycler_module.set_lid_temperature(70)
    thermocycler_module.execute_profile(steps=[
        {'temperature': 4, 'hold_time_minutes': 1},
        {'temperature': 30, 'hold_time_minutes': Enzyme_digestion_time},
        {'temperature': 65, 'hold_time_minutes': 15},
        {'temperature': 4, 'hold_time_minutes': 1}
    ], repetitions=1,block_max_volume=60)
    thermocycler_module.set_block_temperature(4, hold_time_seconds=None)
    thermocycler_module.deactivate_lid()
    thermocycler_module.open_lid()