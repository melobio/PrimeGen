#Code block description: This script block is suitable for the first round of purification steps in a two-step experiment. This script segment defines reagents for a molecular biology protocol on the Opentrons OT-2 robot. It configures four liquids: PCR product, TE buffer (to dissolve DNA), ethanol (to clean magnetic beads), and beads (to bind DNA). Each liquid is given a name, description, and display color for easy identification during the protocol execution.
from opentrons.types import Point
metadata = {
    'protocolName': 'mtb',
    'author': 'MGI-X',
    'apiLevel': '2.14'
}

def run(ctx):
    #Reagent definition
    {round}st_pcr_product = ctx.define_liquid(
        name=f"{round}st_PCR_product",
        description=f"{round}st_pcr_product",
        display_color="#00FF00",
        )
    TE_buffer= ctx.define_liquid(
        name="TE buffer",
        description="TE buffer",#Used to dissolve DNA
        display_color="#00FF00",
        )
    Ethonal = ctx.define_liquid(
        name="Ethonal",
        description="Ethonal",#Used to clean magnetic beads
        display_color="#00FF00",
        )
    Beads = ctx.define_liquid(
        name="Beads",
        description="Beads",#Used to bind DNA
        display_color="#00FF00",
        )