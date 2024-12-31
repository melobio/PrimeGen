from opentrons import protocol_api

metadata = {
    'protocolName': '外显子组酶切_POST-PCR产物纯化_9',
    'author': 'MGIX',
    'description': 'Protocol for DNA fragmentation and recovery',
    'apiLevel': '2.14'
}

#移液函数
def transfer_all(pipette, source_plate, dest_plate, transfer_info,mix_flag = 0,mix_value =[]):
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
        pipette.transfer(volume, temp_source, temp_dest, new_tip='never')
        if mix_flag==1:
            for well in temp_dest:
                pipette.mix(mix_value[0], mix_value[1], well)
        pipette.drop_tip()

def run(protocol: protocol_api.ProtocolContext):

    # Load 200ul filter tip rack
    tips200_1 = protocol.load_labware('opentrons_96_filtertiprack_200ul', '4')
    tips200_2 = protocol.load_labware('opentrons_96_filtertiprack_200ul', '5')
    tips200_3 = protocol.load_labware('opentrons_96_filtertiprack_200ul', '6')
    tips200_4 = protocol.load_labware('opentrons_96_filtertiprack_200ul', '7')
    tips200_5 = protocol.load_labware('opentrons_96_filtertiprack_200ul', '8')
    tips200_6 = protocol.load_labware('opentrons_96_filtertiprack_200ul', '9')
    tips200_7 = protocol.load_labware('opentrons_96_filtertiprack_200ul', '10')
    tips200_8 = protocol.load_labware('opentrons_96_filtertiprack_200ul', '11')

    plate1 = protocol.load_labware('biorad_96_wellplate_200ul_pcr', '3')
    # Load magnetic module and load PCR plate on magnetic module
    mag_module = protocol.load_module('magnetic module', '1')
    mag_plate = mag_module.load_labware('biorad_96_wellplate_200ul_pcr')

    reservoir_plate = protocol.load_labware('usascientific_12_reservoir_22ml', '2')

    # Pipettes
    right_pipette = protocol.load_instrument('p300_multi_gen2', mount='right', tip_racks=[tips200_1,tips200_2,tips200_3,tips200_4,tips200_5,tips200_6,tips200_7,tips200_8])

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
    # Load liquids
    for well in reservoir_plate.columns_by_name()['1']:
        reservoir_plate.wells_by_name()[well.display_name].load_liquid(liquid=DNA_clean_Beads, volume=10400)

    for well in reservoir_plate.columns_by_name()['2']:
        reservoir_plate.wells_by_name()[well.display_name].load_liquid(liquid=TE_Buffer, volume=3328)
    for well in reservoir_plate.columns_by_name()['3']:
        reservoir_plate.wells_by_name()[well.display_name].load_liquid(liquid=Ethanol_80, volume=20800)
    for well in reservoir_plate.columns_by_name()['4']:
        reservoir_plate.wells_by_name()[well.display_name].load_liquid(liquid=Ethanol_80, volume=20800)

    well_name_list = [well.display_name for well in mag_plate.wells()]
    for well_name in well_name_list:
        mag_plate.wells_by_name()[well_name].load_liquid(liquid=PCR_products, volume=100)

    # Execute instructions

    # Transfer
    transfer_times = 12
    transfer_info_right1 = []
    for ii in range(transfer_times):
        transfer_info_right1.append(('1',str(ii + 1), 100))

    transfer_all(right_pipette,reservoir_plate,mag_plate,transfer_info_right1,mix_flag=1,mix_value=[10,50])
    protocol.delay(minutes=5)
    mag_module.engage(height=5)

    transfer_info_right2 = []
    for ii in range(transfer_times):
        transfer_info_right2.append((str(ii + 1),'5',200))

    transfer_all(right_pipette,  mag_plate,reservoir_plate,transfer_info_right2)

    transfer_info_right3 = []
    for ii in range(transfer_times):
        transfer_info_right3.append(('3', str(ii + 1), 200))

    transfer_all(right_pipette, reservoir_plate, mag_plate, transfer_info_right3)
    protocol.delay(seconds=30)

    transfer_info_right4 = []
    for ii in range(transfer_times):
        transfer_info_right4.append((str(ii + 1), '6', 200))

    transfer_all(right_pipette, mag_plate, reservoir_plate, transfer_info_right4)

    transfer_info_right5 = []
    for ii in range(transfer_times):
        transfer_info_right5.append(('4', str(ii + 1), 200))

    transfer_all(right_pipette, reservoir_plate, mag_plate, transfer_info_right5)
    # Pause for 30 seconds
    protocol.delay(seconds=30)
    transfer_info_right6 = []
    for ii in range(transfer_times):
        transfer_info_right6.append((str(ii + 1), '7', 200))

    transfer_all(right_pipette, mag_plate, reservoir_plate, transfer_info_right6)
    protocol.delay(minutes=5)
    # Disengage magnetic module
    mag_module.disengage()
    # Transfer
    transfer_info_right7 = []
    for ii in range(transfer_times):
        transfer_info_right7.append(('2', str(ii + 1), 32))

    transfer_all(right_pipette, reservoir_plate, mag_plate, transfer_info_right7,mix_flag=1,mix_value=[10,20])
    protocol.delay(minutes=5)
    mag_module.engage(height=5)
    protocol.delay(minutes=5)
    # Transfer
    transfer_info_right8 = []
    for ii in range(transfer_times):
        transfer_info_right8.append((str(ii + 1), str(ii + 1), 30))

    transfer_all(right_pipette, mag_plate,plate1, transfer_info_right8)
