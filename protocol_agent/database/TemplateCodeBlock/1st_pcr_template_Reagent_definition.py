#Code block description: This script block is suitable for the first round of PCR reaction steps in a two-step experiment. This script for the Opentrons OT-2 robot defines various liquids to be used in a molecular biology protocol. It configures liquids like DNA template, barcode, primer pools, and PCR enzyme mix, each associated with a unique name, description, and display color. This setup is essential for the automation of loading these reagents in subsequent steps of the protocol.

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
    primer_pools = ctx.define_liquid(
    name="primer pools",
    description="primer pools", 
    display_color="#00FF00", 
    )
    PCR_Enzyme_Mix = ctx.define_liquid(
    name=f"{PCR_round}st_PCR_Enzyme_Mix",
    description=f"{PCR_round}st PCR Enzyme Mix", 
    display_color="#00FF00",
    )
