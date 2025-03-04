from opentrons import protocol_api
## This code block is suitable for magnetic bead double selection steps after enzyme digestion in protocols such as MGIEasy fast/MGIEasy PCRfree etc.
##
metadata = {
    'protocolName': 'MGIEasy fast-2',
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

def run(protocol: protocol_api.ProtocolContext):
    Sample_nums = user_sample_nums
    Breads_volume = input_DNA_Beads_volume
    TE_volume = input_TE_Buffer_volume
    Breads_transfer_volume_1 = input_Breads_first_transfer_volume
    Breads_transfer_volume_2 = input_Breads_second_transfer_volume
    TE_Buffer_transfer_volume = TE_volume/(Sample_nums+8)
    use_col_nums = int(Sample_nums/8)
    if Sample_nums%8 !=0:
        use_col_nums = use_col_nums+1
    # Load 20ul tip rack
    tips20 = protocol.load_labware('opentrons_96_tiprack_20ul', '8')
    # Load 200ul filtered tip racks
    tips200_1 = protocol.load_labware('opentrons_96_tiprack_300ul', '5')
    tips200_2 = protocol.load_labware('opentrons_96_tiprack_300ul', '6')
    tips200_3 = protocol.load_labware('opentrons_96_tiprack_300ul', '7')
    tips200_4 = protocol.load_labware('opentrons_96_tiprack_300ul', '10')
    tips200_5 = protocol.load_labware('opentrons_96_tiprack_300ul', '11')
    deep_plate = protocol.load_labware('nest_96_wellplate_2ml_deep', '2')
    plate1 = protocol.load_labware('biorad_96_wellplate_200ul_pcr', '3')
    plate2 = protocol.load_labware('biorad_96_wellplate_200ul_pcr', '4')
    
    # Load magnetic module and PCR plate on magnetic module
    mag_module = protocol.load_module('magnetic module', '1')
    mag_adapter = mag_module.load_adapter("opentrons_96_flat_bottom_adapter")
    mag_module.open_labware_latch()
    # Pipettes
    left_pipette = protocol.load_instrument('p20_multi_gen2', mount='left', tip_racks=[tips20])
    right_pipette = protocol.load_instrument('p300_multi_gen2', mount='right', tip_racks=[tips200_1,tips200_2,tips200_3,tips200_4,tips200_5])
    left_pipette.flow_rate.aspirate = 3.78  # Set aspirate speed to 50 μL/s
    left_pipette.flow_rate.dispense = 7.56 
    right_pipette.flow_rate.aspirate = 46.43  # Set aspirate speed to 50 μL/s
    right_pipette.flow_rate.dispense = 92.86 
    # Define liquids
    # labeling liquids in wells
    Enzyme_digestion_solution = protocol.define_liquid(
        name="Enzyme digestion reaction solution",
        description="Enzyme digestion reaction solution",
        display_color="#00FF00",
    )
    En_Breads = protocol.define_liquid(
        name="En_Breads",
        description="En_Breads",
        display_color="#00FF00",
    )
    En_TE = protocol.define_liquid(
        name="En_TE",
        description="En_TE",
        display_color="#00FF00",
    )
    # Load liquids
    all_plate_well = [f"{chr(65 + i)}{j}" for j in range(1, 12 + 1) for i in range(8)]
    liquid_list = calculate_liquid(Sample_nums,Breads_volume)
    for ii in range(8):
        well_name = all_plate_well[ii]
        deep_plate.wells_by_name()[well_name].load_liquid(liquid=En_Breads, volume=liquid_list[ii])
        
    liquid_list = calculate_liquid(Sample_nums,TE_volume)
    for ii in range(8):
        well_name = all_plate_well[8:16][ii]
        deep_plate.wells_by_name()[well_name].load_liquid(liquid=En_TE, volume=liquid_list[ii])

    use_plate_well = all_plate_well[:Sample_nums]
    for well_name in use_plate_well:
        plate2.wells_by_name()[well_name].load_liquid(liquid=Enzyme_digestion_solution, volume=60)
    # Execute commands
    # Transfer
    transfer_times = use_col_nums
    transfer_info_right1 = []
    for ii in range(transfer_times):
        transfer_info_right1.append(('1', str(ii + 1), Breads_transfer_volume_1))

    transfer_all(right_pipette, deep_plate, plate2, transfer_info_right1, mix_flag=1, mix_value=[10, 50])

    protocol.delay(minutes=5)
    protocol.move_labware(labware=plate2, new_location=mag_adapter)
    mag_module.engage(height_from_base=5)
    protocol.delay(minutes=5)

    transfer_info_right2 = []
    for ii in range(transfer_times):
        transfer_info_right2.append((str(ii + 1), str(ii + 1), 60+Breads_transfer_volume_1))

    transfer_all(right_pipette, plate2, plate1, transfer_info_right2)
    mag_module.disengage()

    transfer_info_left1 = []
    for ii in range(transfer_times):
        transfer_info_left1.append(('1', str(ii + 1), Breads_transfer_volume_2))

    transfer_all(left_pipette, deep_plate, plate1, transfer_info_left1)
    # mix
    for ii in range(transfer_times):
        dest = plate1.columns_by_name()[str(ii + 1)]
        right_pipette.pick_up_tip()
        right_pipette.mix(10, 50, dest)
        right_pipette.drop_tip()

    protocol.delay(minutes=5)
    protocol.move_labware(labware=plate2, new_location='4')
    protocol.move_labware(labware=plate1, new_location=mag_adapter)
    mag_module.engage(height_from_base=5)
    protocol.delay(minutes=5)
    # Transfer
    transfer_info_right3 = []
    for ii in range(transfer_times):
        transfer_info_right3.append((str(ii + 1), '3', 108))

    transfer_all(right_pipette, plate1, deep_plate, transfer_info_right3)
    mag_module.disengage()
    protocol.move_labware(labware=plate1, new_location='3')

    transfer_info_right4 = []
    for ii in range(transfer_times):
        transfer_info_right4.append(('2', str(ii + 1), TE_Buffer_transfer_volume))

    transfer_all(right_pipette, deep_plate, plate1, transfer_info_right4, mix_flag=1, mix_value=[10, 50])

