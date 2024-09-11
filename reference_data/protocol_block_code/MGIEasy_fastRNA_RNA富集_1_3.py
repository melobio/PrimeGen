from opentrons import protocol_api

metadata = {
    'protocolName': 'MGIEasy_fastRNA_RNA富集_1_3',
    'author': 'MGIX',
    'description': 'Protocol for DNA fragmentation and recovery',
    'apiLevel': '2.14'
}


# 移液函数
def transfer_all(pipette, source_plate, dest_plate, transfer_info, mix_flag=0, mix_value=[]):
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
    Sample_nums = user_sample_nums
    TRIS_HCL_buffer_volume = input_TRIS_HCL_buffer_volume_2
    use_col_nums = int(Sample_nums / 8)
    if Sample_nums % 8 != 0:
        use_col_nums = use_col_nums + 1
    # 加载200ul吸头架
    tips20_1 = protocol.load_labware('opentrons_96_tiprack_20ul', '7')
    tips20_2 = protocol.load_labware('opentrons_96_tiprack_20ul', '8')
    tips200_1 = protocol.load_labware('opentrons_96_tiprack_300ul', '5')

    temp_shaker = protocol.load_module('temperature module gen2', '1')
    plate1 = protocol.load_labware('biorad_96_wellplate_200ul_pcr', '6')
    # 加载磁性模块和在磁性模块上加载PCR板
    mag_module = protocol.load_module('magnetic module', '3')
    mag_plate = mag_module.load_labware('biorad_96_wellplate_200ul_pcr')
    deep_plate = protocol.load_labware('nest_96_wellplate_2ml_deep', '9')
    # Pipettes
    left_pipette = protocol.load_instrument('p20_multi_gen2', mount='left', tip_racks=[tips20_1, tips20_2])
    right_pipette = protocol.load_instrument('p300_multi_gen2', mount='right', tip_racks=[tips200_1])
    right_pipette.flow_rate.aspirate = 46.43  # Set aspirate speed to 50 μL/s
    right_pipette.flow_rate.dispense = 92.86
    left_pipette.flow_rate.aspirate = 3.78  # Set aspirate speed to 50 μL/s
    left_pipette.flow_rate.dispense = 7.56
    # 定义液体
    # labeling liquids in wells
    RNA_enrichment_bead_washbuffer = protocol.define_liquid(
        name="RNA_enrichment_bead_washbuffer",
        description="RNA_enrichment_bead_washbuffer",
        display_color="#00FF00",
    )
    TRIS_HCL = protocol.define_liquid(
        name="TRIS-HCL",
        description="TRIS-HCL",
        display_color="#00FF00",
    )

    # 加载液体
    all_plate_well = [f"{chr(65 + i)}{j}" for j in range(1, 12 + 1) for i in range(8)]
    liquid_list_1 = calculate_liquid(Sample_nums, TRIS_HCL_buffer_volume)
    for ii in range(8):
        deep_plate.wells_by_name()[all_plate_well[ii]].load_liquid(liquid=TRIS_HCL, volume=liquid_list_1[ii])

    use_plate_well = all_plate_well[:Sample_nums]
    for well_name in use_plate_well:
        mag_plate.wells_by_name()[well_name].load_liquid(liquid=RNA_enrichment_bead_washbuffer, volume=50)

        # 执行指令
    mag_module.engage(height_from_base=5)
    protocol.delay(minutes=2)

    # 移液
    transfer_times = use_col_nums
    transfer_info_right1 = []
    for ii in range(transfer_times):
        transfer_info_right1.append((str(ii + 1),'2', 50))

    transfer_all(right_pipette, mag_plate,deep_plate,  transfer_info_right1)
    mag_module.disengage()

    transfer_info_left1 = []
    for ii in range(transfer_times):
        transfer_info_left1.append(('1',str(ii + 1), 12))
    transfer_all(left_pipette,deep_plate, mag_plate,transfer_info_left1,mix_flag=1,mix_value=[10,5])

    temp_shaker.deactivate()
    temp_shaker.open_labware_latch()
    protocol.move_labware(labware=mag_plate, new_location='1')
    temp_shaker.set_temperature(80)

    temp_shaker.close_labware_latch()
    protocol.delay(minutes=2)
    temp_shaker.deactivate()
    temp_shaker.open_labware_latch()
    protocol.move_labware(labware=mag_plate, new_location='3')
    mag_module.engage(height_from_base=5)
    protocol.delay(minutes=2)
    transfer_info_left2 = []
    for ii in range(transfer_times):
        transfer_info_left2.append((str(ii + 1),str(ii + 1),  10))
    transfer_all(left_pipette, mag_plate,plate1, transfer_info_left2)


