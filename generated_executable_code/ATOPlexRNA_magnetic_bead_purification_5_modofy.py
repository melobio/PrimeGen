from opentrons import protocol_api

metadata = {
    'protocolName': 'ATOPlexDNA-4',
    'author': 'MGIX',
    'description': 'Protocol for DNA fragmentation and recovery',
    'apiLevel': '2.14'
}


#移液函数
def transfer_all(pipette, source_plate, dest_plate, transfer_info,mix_flag = 0,mix_value =[]):
    """
    General function for transferring samples.
    Args:
    - pipette: the pipette to use
    - source_plate: source plate
    - dest_plate: destination plate
    - transfer_info: a list containing transfer information, each element is a tuple (source column name, destination column name, volume)
    - mix_value: mixing information, a list, [number, volume]
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
def get_beads_volume(Barcode_type,panel_amp_length):
    beads_volume_1st = 0
    beads_volume_2nd = 0
    if Barcode_type == 'single_pcr_barcode':
        if panel_amp_length >=100 and panel_amp_length <200:
            beads_volume_1st = 32.5
            beads_volume_2nd = 27.5
        elif panel_amp_length >=200 and panel_amp_length <300:
            beads_volume_1st = 25
            beads_volume_2nd = 22.5
        else:
            beads_volume_1st = 22.5
            beads_volume_2nd = 20
    elif Barcode_type == 'double_pcr_barcode':
        if panel_amp_length >=60 and panel_amp_length <100:
            beads_volume_1st = 37.5
            beads_volume_2nd = 37.5
        elif panel_amp_length >=100 and panel_amp_length <200:
            beads_volume_1st = 32.5
            beads_volume_2nd = 25
        elif panel_amp_length >=200 and panel_amp_length <300:
            beads_volume_1st = 25
            beads_volume_2nd = 20
        else:
            beads_volume_1st = 22.5
            beads_volume_2nd = 20

    return beads_volume_1st,beads_volume_2nd


def run(protocol: protocol_api.ProtocolContext):
    Sample_nums = 60
    Beads_volume = 1360
    TE_Buffer_volume = 1700
    Ethanol_80_volume = 21760
    use_col_nums = int(Sample_nums/8)
    if Sample_nums%8 !=0:
        use_col_nums = use_col_nums+1

    Spike_in_Control_PCR_block = false
    Barcode_type = 'single_pcr_barcode'
    panel_amp_length = 300
    
    beads_volume_1st,beads_volume_2nd = get_beads_volume(Barcode_type,panel_amp_length)
    
    if Spike_in_Control_PCR_block:
        beads_volume_2nd = beads_volume_2nd +1

    tips200_1 = protocol.load_labware('opentrons_96_tiprack_300ul', '4')
    tips200_2 = protocol.load_labware('opentrons_96_tiprack_300ul', '5')
    tips200_3 = protocol.load_labware('opentrons_96_tiprack_300ul', '6')
    tips200_4 = protocol.load_labware('opentrons_96_tiprack_300ul', '7')
    tips200_5 = protocol.load_labware('opentrons_96_tiprack_300ul', '8')
    tips200_6 = protocol.load_labware('opentrons_96_tiprack_300ul', '9')
    tips200_7 = protocol.load_labware('opentrons_96_tiprack_300ul', '10')
    tips200_8 = protocol.load_labware('opentrons_96_tiprack_300ul', '11')


    plate1 = protocol.load_labware('biorad_96_wellplate_200ul_pcr', '3')
    #
    mag_module = protocol.load_module('magnetic module gen2', '1')
    mag_plate = mag_module.load_labware('biorad_96_wellplate_200ul_pcr')

    reservoir_plate = protocol.load_labware('usascientific_12_reservoir_22ml', '2')

    # Pipettes
    right_pipette = protocol.load_instrument('p300_multi_gen2', mount='right', tip_racks=[tips200_1,tips200_2,tips200_3,tips200_4,tips200_5,tips200_6,tips200_7,tips200_8])
    right_pipette.flow_rate.aspirate = 46.43  # Set aspirate speed to 50 μL/s
    right_pipette.flow_rate.dispense = 92.86 
    #
    # labeling liquids in wells
    PCR_2nd = protocol.define_liquid(
        name="PCR_2nd",
        description="PCR_2nd",
        display_color="#00FF00",
    )
    DNA_clean_Beads = protocol.define_liquid(
        name="DNA_clean_Beads",
        description="DNA_clean_Beads",
        display_color="#00FF00",
    )
    TE_Buffer = protocol.define_liquid(
        name="TE_Buffer",
        description="TE_Buffer",
        display_color="#00FF00",
    )
    Ethanol_80 = protocol.define_liquid(
        name="Ethanol_80",
        description="Ethanol_80",
        display_color="#00FF00",
    )
    #
    all_plate_well = [f"{chr(65 + i)}{j}" for j in range(1, 12 + 1) for i in range(8)]
    reservoir_plate.wells_by_name()['A1'].load_liquid(liquid=DNA_clean_Beads, volume=Beads_volume)
    reservoir_plate.wells_by_name()['A2'].load_liquid(liquid=TE_Buffer, volume=TE_Buffer_volume)
    reservoir_plate.wells_by_name()['A3'].load_liquid(liquid=Ethanol_80, volume=Ethanol_80_volume/2)
    reservoir_plate.wells_by_name()['A4'].load_liquid(liquid=Ethanol_80, volume=Ethanol_80_volume/2)

    use_plate_well = all_plate_well[:Sample_nums]
    for well_name in use_plate_well:
        mag_plate.wells_by_name()[well_name].load_liquid(liquid=PCR_2nd, volume=beads_volume_2nd)

    transfer_times = use_col_nums
    transfer_info_right1 = []
    for ii in range(transfer_times):
        transfer_info_right1.append(('1', str(ii + 1), beads_volume_2nd))

    transfer_all(right_pipette, reservoir_plate, mag_plate, transfer_info_right1, mix_flag=1, mix_value=[10, 25])
    protocol.delay(minutes=5)
    mag_module.engage(height_from_base=5)
    protocol.delay(minutes=5)

    transfer_info_right2 = []
    for ii in range(transfer_times):
        transfer_info_right2.append((str(ii + 1), '5', 2*beads_volume_2nd))

    transfer_all(right_pipette, mag_plate, reservoir_plate, transfer_info_right2)

    transfer_info_right3 = []
    for ii in range(transfer_times):
        transfer_info_right3.append(('3', str(ii + 1), 160))

    transfer_all(right_pipette, reservoir_plate, mag_plate, transfer_info_right3)
    protocol.delay(seconds=30)

    transfer_info_right4 = []
    for ii in range(transfer_times):
        transfer_info_right4.append((str(ii + 1), '6', 160))

    transfer_all(right_pipette, mag_plate, reservoir_plate, transfer_info_right4)

    transfer_info_right5 = []
    for ii in range(transfer_times):
        transfer_info_right5.append(('4', str(ii + 1), 160))

    transfer_all(right_pipette, reservoir_plate, mag_plate, transfer_info_right5)

    protocol.delay(seconds=30)

    transfer_info_right6 = []
    for ii in range(transfer_times):
        transfer_info_right6.append((str(ii + 1), '7', 160))

    transfer_all(right_pipette, mag_plate, reservoir_plate, transfer_info_right6)
    protocol.delay(minutes=5)

    mag_module.disengage()

    transfer_info_right7 = []
    for ii in range(transfer_times):
        transfer_info_right7.append(('2', str(ii + 1), 25))

    transfer_all(right_pipette, reservoir_plate, mag_plate, transfer_info_right7)

    protocol.delay(minutes=5)
    mag_module.engage(height_from_base=5)
    protocol.delay(minutes=5)

    transfer_info_right8= []
    for ii in range(transfer_times):
        transfer_info_right8.append((str(ii + 1), str(ii + 1), 23))

    transfer_all(left_pipette, mag_plate, plate1, transfer_info_right8)
