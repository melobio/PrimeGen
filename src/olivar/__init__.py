#!/usr/bin/env python3
# -*- coding: utf-8 -*-


'''
Olivar is a Python3 software for multiplex PCR tiling design. 
Olivar implements a novel algorithm to reduce non-specific amplifications in PCR, 
while avoiding primer design at SNPs and other undesired regions at the same time. 
Olivar also optimize for primer dimers with the SADDLE algorithm (https://www.nature.com/articles/s41467-022-29500-4). 
'''


__author__ = 'Michael X. Wang'
__version__ = '1.0.2'


import os
import sys
# CURRENT_DIR = os.path.dirname(__file__)
# print(CURRENT_DIR)
# print(sys.path)
# sys.path.append(CURRENT_DIR)

#from pipeline import build, tiling, visualize, validate
#from .pipeline_snp import build, tiling, visualize, validate
#from .pipeline_ada import build, tiling, visualize, validate
# from .pipeline_snp_species import build, tiling, visualize, validate
# from .pipeline_snp_species_round2_1tube_primer3 import build, tiling, visualize, validate
from .pipeline_snp_species_round2_1tube_primer3_get_primer import build, tiling, visualize, validate
