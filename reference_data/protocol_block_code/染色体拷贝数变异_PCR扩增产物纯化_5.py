from opentrons import protocol_api

metadata = {
    'protocolName': 'Chromosome Copy Number Variation PCR Product Purification 5',
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
    - mix_value: Mixing information as a list [number of times, volume]
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


def run(protocol: protocol_api.ProtocolContext):
    Sample_nums = user_sample_nums
    Beads_volume =input_DNA_Beads_volume
    use_col_nums = int(Sample_nums / 8)
    if Sample_nums % 8 != 0:
        use_col_nums = use_col_nums + 1

    # Load 20ul tip rack
    tips20 = protocol.load_labware('opentrons_96_tiprack_20ul', '6')
    # Load 200ul filtered tip racks
    tips200_1 = protocol.load_labware('opentrons_96_tiprack_300ul', '4')
    tips200_2 = protocol.load_labware('opentrons_96_tiprack_300ul', '5')
    tips200_3 = protocol.load_labware('opentrons_96_tiprack_300ul', '7')

    plate1 = protocol.load_labware('biorad_96_wellplate_200ul_pcr', '3')
    # Load magnetic module and PCR plate on magnetic module
    mag_module = protocol.load_module('magnetic module gen2', '1')
    mag_plate = mag_module.load_labware('biorad_96_wellplate_200ul_pcr')

    reservoir_plate = protocol.load_labware('usascientific_12_reservoir_22ml', '2')

    # Pipettes
    left_pipette = protocol.load_instrument('p20_multi_gen2', mount='left', tip_racks=[tips20])
    right_pipette = protocol.load_instrument('p300_multi_gen2', mount='right',
                                             tip_racks=[tips200_1, tips200_2, tips200_3])
    left_pipette.flow_rate.aspirate = 3.78  # Set aspirate speed to 50 μL/s
    left_pipette.flow_rate.dispense = 7.56
    right_pipette.flow_rate.aspirate = 46.43  # Set aspirate speed to 50 μL/s
    right_pipette.flow_rate.dispense = 92.86
    # Define liquids
    # labeling liquids in wells
    PCR_products = protocol.define_liquid(
        name="PCR_products",
        description="PCR_products",
        display_color="#00FF00",
    )
    DNA_clean_Beads = protocol.define_liquid(
        name="DNA_clean_Beads",
        description="DNA_clean_Beads",
        display_color="#00FF00",
    )
    # Load liquids
    all_plate_well = [f"{chr(65 + i)}{j}" for j in range(1, 12 + 1) for i in range(8)]
    reservoir_plate.wells_by_name()['A1'].load_liquid(liquid=DNA_clean_Beads, volume=Beads_volume)

    use_plate_well = all_plate_well[:Sample_nums]
    for well_name in use_plate_well:
        mag_plate.wells_by_name()[well_name].load_liquid(liquid=PCR_products, volume=50)

    # Execute commands
    # Transfer
    transfer_times = use_col_nums
    transfer_info_right1 = []
    for ii in range(transfer_times):
        transfer_info_right1.append(('1', str(ii + 1), 35))

    transfer_all(right_pipette, reservoir_plate, mag_plate, transfer_info_right1, mix_flag=1, mix_value=[10, 50])
    # Pause for 5 minutes
    protocol.delay(minutes=5)
    mag_module.engage(height_from_base=5)
    protocol.delay(minutes=5)

    transfer_info_right2 = []
    for ii in range(transfer_times):
        transfer_info_right2.append((str(ii + 1), str(ii + 1), 85))

    transfer_all(right_pipette, mag_plate, plate1, transfer_info_right2)

    transfer_info_left1 = []
    for ii in range(transfer_times):
        transfer_info_left1.append(('1', str(ii + 1), 15))

    transfer_all(left_pipette, reservoir_plate, plate1, transfer_info_left1)

    for ii in range(transfer_times):
        temp_dest = plate1.columns_by_name()[str(ii + 1)]
        right_pipette.pick_up_tip()
        right_pipette.mix(10, 50, temp_dest)
        right_pipette.drop_tip()
    protocol.delay(minutes=5)