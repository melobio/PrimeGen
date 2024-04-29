#Code block description: This script block is suitable for the second round of purification steps in a two-step experiment. This script defines various reagents for use in a molecular biology protocol on an Opentrons OT-2 robot, including PCR products, TE buffer, ethanol, and beads. Each liquid is specified with a name, purpose, and display color to facilitate accurate and automated handling during experimental procedures.

from opentrons.types import Point
metadata = {
    'protocolName': 'mtb',
    'author': 'MGI-X',
    'apiLevel': '2.14'
}

def run(ctx):
    #Reagent definition
    {round}nd_pcr_product = ctx.define_liquid(
        name=f"{round}nd_PCR_product",
        description=f"{round}nd_pcr_product",
        display_color="#00FF00",
        )
    {round}nd_TE_buffer= ctx.define_liquid(
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