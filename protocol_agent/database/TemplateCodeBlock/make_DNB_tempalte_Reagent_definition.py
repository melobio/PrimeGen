#Code block description: This script block is suitable for the make DNB steps in a two-step experiment. This script defines various reagents for a DNA library preparation protocol on the Opentrons OT-2 robot. It includes Low TE buffer, dsDNA library, DNB buffer, two enzyme mixes for DNB formation, and a stop buffer. Each reagent is characterized by name, description, and display color for automated handling.

from opentrons.types import Point
metadata = {
    'protocolName': 'mtb',
    'author': 'MGI-X',
    'apiLevel': '2.14'
}

def run(ctx):
    #Reagent definition
    Low_TE_buffer = ctx.define_liquid(
        name="Low_TE_buffer",
        description="Low_TE_buffer",#Low_TE_buffer
        display_color="#00FF00",
        )
    dsDNA_libarary = ctx.define_liquid(
        name="dsDNA_libarary",
        description="dsDNA_libarary",#dsDNA_libarary
        display_color="#00FF00",
        )
    make_DNB_buffer = ctx.define_liquid(
        name="make_DNB_buffer",
        description="make_DNB_buffer",#make DNB buffer
        display_color="#00FF00",
        )
    make_DNB_enzyme_mix_1= ctx.define_liquid(
        name="make_DNB_enzyme_mix_1",
        description="make_DNB_enzyme_mix_1",#make_DNB_enzyme_mix_1
        display_color="#00FF00",
        )
    make_DNB_enzyme_mix_2= ctx.define_liquid(
        name="make_DNB_enzyme_mix_2",
        description="make_DNB_enzyme_mix_2",#make_DNB_enzyme_mix_2
        display_color="#00FF00",
        )
    make_DNB_stop_buffer = ctx.define_liquid(
        name="Beads",
        description="Beads",#Terminate DNB reaction
        display_color="#00FF00",
        )