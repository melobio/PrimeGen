#Code block description: This script block is suitable for the second round of PCR reaction steps in a two-step experiment. This script defines various liquids for a PCR experiment on an Opentrons OT-2 robot, including DNA template, barcode, UR, and PCR enzyme mix. Each liquid is given a name, description, and display color to facilitate automated handling and tracking during the protocol execution.

from opentrons.types import Point
metadata = {
    'protocolName': 'mtb',
    'author': 'MGI-X',
    'apiLevel': '2.14'
}

def run(ctx):
    #Reagent definition
    DNA_TEMPLATE = ctx.define_liquid(
    name=f"DNA_TEMPLATE_{PCR_round}st_PCR",
    description=f"DNA_TEMPLATE_{PCR_round}st_PCR", 
    display_color="#00FF00",
    )
    Barcode_{PCR_round} = ctx.define_liquid(
    name=f"Barcode_{PCR_round}",
    description="Barcode_{PCR_round}", 
    display_color="#00FF00",
    )
    UR = ctx.define_liquid(
    name="UR",
    description="UR",
    display_color="#00FF00",
    )
    PCR_Enzyme_Mix = ctx.define_liquid(
    name=f"{PCR_round}st_PCR_Enzyme_Mix",
    description=f"{PCR_round}st PCR Enzyme Mix", 
    display_color="#00FF00",
    )