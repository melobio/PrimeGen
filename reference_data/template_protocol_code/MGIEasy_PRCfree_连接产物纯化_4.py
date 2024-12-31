from opentrons import protocol_api

metadata = {
    'protocolName': 'MGIEasy_PRCfree_连接产物纯化_5',
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


def run(protocol: protocol_api.ProtocolContext):
    Sample_nums = input_sample_nums_1
    Beads_volume = input_DNA_Beads_volume_4
    TE_Buffer_volume = input_TE_Buffer_volume_4
    use_col_nums = int(Sample_nums / 8)
    if Sample_nums % 8 != 0:
        use_col_nums = use_col_nums + 1
    # Load 20ul tip rack
    #tips20 = protocol.load_labware('opentrons_96_tiprack_20ul', '9')
    # Load 200ul filter tip rack
    tips200_1 = protocol.load_labware('opentrons_96_tiprack_300ul', '4')
    tips200_2 = protocol.load_labware('opentrons_96_tiprack_300ul', '5')
    tips200_3 = protocol.load_labware('opentrons_96_tiprack_300ul', '6')
    tips200_4 = protocol.load_labware('opentrons_96_tiprack_300ul', '7')
    tips200_5 = protocol.load_labware('opentrons_96_tiprack_300ul', '8')

    plate1 = protocol.load_labware('biorad_96_wellplate_200ul_pcr', '3')
    # Load magnetic module and load PCR plate on magnetic module
    mag_module = protocol.load_module('magnetic module gen2', '1')
    mag_adapter = mag_module.load_adapter("opentrons_96_flat_bottom_adapter")
    mag_module.open_labware_latch()

    reservoir_plate = protocol.load_labware('usascientific_12_reservoir_22ml', '2')
    # Pipettes
    right_pipette = protocol.load_instrument('p300_multi_gen2', mount='right', tip_racks=[tips200_1, tips200_2, tips200_3, tips200_4, tips200_5])
    right_pipette.flow_rate.aspirate = 46.43  # Set aspirate speed to 50 μL/s
    right_pipette.flow_rate.dispense = 92.86 
    # Define liquids
    connect_products = protocol.define_liquid(
        name="Connect products",
        description="Connect products",
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
    reservoir_plate["A1"].load_liquid(liquid=En_Breads, volume=TE_Buffer_volume)
    reservoir_plate["A2"].load_liquid(liquid=En_TE, volume=Beads_volume)
    use_plate_well = all_plate_well[:Sample_nums]
    for well_name in use_plate_well:
        plate1.wells_by_name()[well_name].load_liquid(liquid=connect_products, volume=80)

    # Execute instructions
    # Transfer
    transfer_times = use_col_nums
    transfer_info_right1 = []
    for ii in range(transfer_times):
        transfer_info_right1.append(('1', str(ii + 1), 20))
    transfer_all(right_pipette, reservoir_plate, plate1, transfer_info_right1, mix_flag=1, mix_value=[10, 50])

    transfer_info_right2 = []
    for ii in range(transfer_times):
        transfer_info_right2.append(('2', str(ii + 1), 20))
    transfer_all(right_pipette, reservoir_plate, plate1, transfer_info_right2, mix_flag=1, mix_value=[10, 50])
    # Pause for 5 minutes
    protocol.delay(minutes=5)
    protocol.move_labware(labware=plate1, new_location=mag_adapter)
    mag_module.engage(height_from_base=5)
    protocol.delay(minutes=5)

    transfer_info_right3 = []
    for ii in range(transfer_times):
        transfer_info_right3.append((str(ii + 1), '5', 120))
    transfer_all(right_pipette, plate1, reservoir_plate, transfer_info_right3)
    protocol.move_labware(labware=plate1, new_location='4')

    transfer_info_right4 = []
    for ii in range(transfer_times):
        transfer_info_right4.append(('1', str(ii + 1), 30))
    transfer_all(right_pipette, reservoir_plate, plate1, transfer_info_right4)

    transfer_info_right5 = []
    for ii in range(transfer_times):
        transfer_info_right5.append(('2', str(ii + 1), 20))

    transfer_all(right_pipette, reservoir_plate, plate1, transfer_info_right5, mix_flag=1, mix_value=[10, 25])
    protocol.delay(minutes=5)
