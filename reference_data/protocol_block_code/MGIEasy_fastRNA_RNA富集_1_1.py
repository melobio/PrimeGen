from opentrons import protocol_api

metadata = {
    'protocolName': 'MGIEasy_fastRNA_RNA富集_1_1',
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
    Washing_Buffer_volume = input_Washing_Buffer_volume_1
    Beads_volume = input_RNA_Enrichment_beads_volume
    use_col_nums = int(Sample_nums / 8)
    if Sample_nums % 8 != 0:
        use_col_nums = use_col_nums + 1

    # Load 200ul tip rack
    tips200_1 = protocol.load_labware('opentrons_96_tiprack_300ul', '5')
    tips200_2 = protocol.load_labware('opentrons_96_tiprack_300ul', '7')
    tips200_3 = protocol.load_labware('opentrons_96_tiprack_300ul', '8')
    tips200_4 = protocol.load_labware('opentrons_96_tiprack_300ul', '10')
    tips200_5 = protocol.load_labware('opentrons_96_tiprack_300ul', '11')

    temp_shaker = protocol.load_module('temperature module gen2', '1')
    temp_plate = temp_shaker.load_labware('biorad_96_wellplate_200ul_pcr')
    # Load thermocycler module, make sure to specify the correct position number
    plate1 = protocol.load_labware('biorad_96_wellplate_200ul_pcr', '6')
    # Load magnetic module and PCR plate on magnetic module
    mag_module = protocol.load_module('magnetic module', '3')
    #mag_plate = mag_module.load_labware('biorad_96_wellplate_200ul_pcr')
    deep_plate = protocol.load_labware('nest_96_wellplate_2ml_deep', '9')
    # Pipettes

    right_pipette = protocol.load_instrument('p300_multi_gen2', mount='right', tip_racks=[tips200_1,tips200_2,tips200_3,tips200_4,tips200_5])
    right_pipette.flow_rate.aspirate = 46.43  # Set aspirate speed to 50 μL/s
    right_pipette.flow_rate.dispense = 92.86
    # Define liquids
    # labeling liquids in wells
    RNA_sample = protocol.define_liquid(
        name="RNA_sample",
        description="RNA_sample",
        display_color="#00FF00",
    )
    Beads = protocol.define_liquid(
        name="Beads",
        description="Beads",
        display_color="#00FF00",
    )
    Washing_Buffer = protocol.define_liquid(
        name="Washing_Buffer",
        description="Washing_Buffer",
        display_color="#00FF00",
    )

    # Load liquids
    all_plate_well = [f"{chr(65 + i)}{j}" for j in range(1, 12 + 1) for i in range(8)]
    liquid_list_1 = calculate_liquid(Sample_nums, Beads_volume)
    liquid_list_2 = calculate_liquid(Sample_nums, Washing_Buffer_volume)

    for ii in range(8):
        deep_plate.wells_by_name()[all_plate_well[ii]].load_liquid(liquid=Beads, volume=liquid_list_1[ii])
        deep_plate.wells_by_name()[all_plate_well[8:16][ii]].load_liquid(liquid=Washing_Buffer, volume=liquid_list_2[ii])


    use_plate_well = all_plate_well[:Sample_nums]
    for well_name in use_plate_well:
        plate1.wells_by_name()[well_name].load_liquid(liquid=RNA_sample, volume=25)

        # Execute instructions
    # Set temperature and shaking parameters
    temp_shaker.set_temperature(65)  # Set temperature to 65°C
    temp_shaker.deactivate()
    temp_shaker.close_labware_latch()
    temp_shaker.await_temperature(65)

    # Liquid transfer
    transfer_times = use_col_nums
    transfer_info_right1 = []
    for ii in range(transfer_times):
        transfer_info_right1.append((str(ii + 1), str(ii + 1), 25))

    transfer_all(right_pipette, plate1,temp_plate,  transfer_info_right1)
    protocol.delay(minutes=5)
    transfer_info_right2 = []
    for ii in range(transfer_times):
        transfer_info_right2.append(('1', str(ii + 1), 25))

    transfer_all(right_pipette, deep_plate,temp_plate,  transfer_info_right2,mix_flag=1,mix_value=[10,30])
    protocol.delay(minutes=5)
    temp_shaker.deactivate()
    temp_shaker.open_labware_latch()

    protocol.move_labware(labware=temp_plate, new_location='3')
    mag_module.engage(height_from_base=5)
    protocol.delay(minutes=2)

    transfer_info_right3 = []
    for ii in range(transfer_times):
        transfer_info_right3.append((str(ii + 1),'5', 50))

    transfer_all(right_pipette,temp_plate,deep_plate, transfer_info_right3)
    mag_module.disengage()

    transfer_info_right4 = []
    for ii in range(transfer_times):
        transfer_info_right4.append(('2',str(ii + 1),  50))
    transfer_all(right_pipette, deep_plate,temp_plate, transfer_info_right4,mix_flag=1,mix_value=[10,30])
    mag_module.engage(height_from_base=5)
    protocol.delay(minutes=2)

    transfer_info_right5 = []
    for ii in range(transfer_times):
        transfer_info_right5.append((str(ii + 1),'5', 50))
    transfer_all(right_pipette,  temp_plate,deep_plate,  transfer_info_right5)

