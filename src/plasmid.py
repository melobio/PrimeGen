# This is a sample Python script.

# Press ⌃R to execute it or replace it with your code.
# Press Double ⇧ to search everywhere for classes, files, tool windows, actions, and settings.

from olivar import build, tiling#, save, validate
import os
import re
import pandas as pd
import ast


#print(os.path.abspath('./'))
candidates_dict = {}

# second round
round2_exclude_dict = {}



if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Primer Design Code')
    parser.add_argument('--fasta_path', help='the path of fasta file')
    parser.add_argument('--out_path', help='the path of output')
    parser.add_argument('--gene_name', help='the gene name')
    parser.add_argument('--max_amp_len', help='Maximum amplicon length', type=int, default=270)
    parser.add_argument('--min_amp_len', help='Minimum amplicon length', type=int, default=210)
    parser.add_argument('--temperature', help='PCR annealing temperature', type=int, default=53)
    parser.add_argument('--min_GC', help='Minimum GC content of a primer', type=float, default=0.4)
    parser.add_argument('--max_GC', help='Maximum GC content of a primer', type=float, default=0.65)
    parser.add_argument('--max_len', help='Maximum length of a primer', type=int, default=23)
    parser.add_argument('--start_pos', help='the start_pos', type=int, default=0)
    parser.add_argument('--end_pos', help='the end_pos', type=int, default=0)
    parser.add_argument('--snp_path', help='the path of snp file', default=None)
    args = parser.parse_args()
    fasta_path = args.fasta_path
    out_path = args.out_path
    var_path = args.snp_path
    start_pos = args.start_pos
    end_pos = args.end_pos
    temperature = args.temperature
    max_amp_len = args.max_amp_len
    min_amp_len = args.min_amp_len
    min_GC = args.min_GC
    max_GC = args.max_GC
    gene_name = args.gene_name
    max_len = args.max_len
    print('max_amp_len', max_amp_len)
    print('min_amp_len', min_amp_len)
    print('min_GC', min_GC)
    print('max_GC', max_GC)
    print('max_len', max_len)
    print('Start primer design')
    build(
        fasta_path = fasta_path,
        var_path = var_path,
        BLAST_db = None,
        ## Must write the database name output from makeblastdb
        #BLAST_db = '/Users/wangyi/PycharmProjects/pcr_1/twt_设计/dt1/dt1_db/dt1_db', ## Must write the database name output from makeblastdb
        #BLAST_db= '/Users/wangyi/PycharmProjects/pcr_1/primer_agent/snp_target_typing/db_NC_000962.3/NC_000962',
        #Optional, path to the BLAST database.
        # Note that this path should end with the name of the BLAST database (e.g., "example_input/Human/GRCh38_primary").
        out_path = out_path,
        title = f'primer_design_{gene_name}', # Name of the Olivar reference file [olivar-ref].
        threads = 1 ,# Number of threads.
        start_pos=start_pos,
        end_pos=end_pos,
    )

    tiling(
        #ref_path='./hsp65/hsp65-ref.olvr',  # Path to the Olivar reference file (.olvr).
        #ref_path='./primer_agent/snp_target_typing/20240119_target_typing_ref.olvr',
        #ref_path='./primer_agent/plasmid/20240123_plasmid_ref.olvr',
        ref_path=f'/reference_data/primer_result/primer_design_{gene_name}.olvr',

        #out_path='./hsp65/design',  # Output directory.
        out_path=out_path,
        #out_path='./twt_设计/dt1/round2/',

        # title='dt1-p-design',  # Name of this design [olivar-design].
        #title='mtb23_no_species-p-design',
        title=f'primer_design_{gene_name}_result',
#         title='dt1_round2_design10',

        # max_amp_len=int(args.max_amp_len),  # Maximum amplicon length [420].
        # min_amp_len=int(args.min_amp_len),  # Minimum amplicon length. 0.9*{max-amp-len} if set as None. Minimum 120.
        max_amp_len=max_amp_len,  # Maximum amplicon length [420].
        min_amp_len=min_amp_len,  # Minimum amplicon length. 0.9*{max-amp-len} if set
        # max_amp_len=270,  # Maximum amplicon length [420].
        # min_amp_len=210,  # Minimum amplicon length. 0.9*{max-amp-len} if set
        w_egc=1,  # Weight for extreme GC content [1.0].
        w_lc=1,  # Weight for low sequence complexity [1.0].
        w_ns=1,  # Weight for non-specificity [1.0]. dt7k:5;
        w_var=1,  # Weight for variations [1.0].

        temperature=temperature,  # PCR annealing temperature [60.0].
        salinity=0.18,  # Concentration of monovalent ions in units of molar [0.18].
        dG_max=-11.8,  # Maximum free energy change of a primer in kcal/mol [-11.8].
        min_GC=min_GC,  # Minimum GC content of a primer [0.2].
        max_GC=max_GC,  # Maximum GC content of a primer [0.75].
        min_complexity=0.1,  # Minimum sequence complexity of a primer [0.4].
        max_len=max_len,  # Maximum length of a primer [36].
        #check_var=True,  # Filter out primer candidates with variations within 5nt of 3' end [False].
        # Setting check_var=True is not recommended when a lot of variations are provided,
        # since this would significantly reduce the number of primer candidates.

        fP_prefix='',  # Prefix of forward primer. Empty string '' by default.
        rP_prefix='',  # Prefix of reverse primer. Empty string '' by default.

        seed=10,  # Random seed for optimizing primer design regions and primer dimer [10].
        threads=10,  # Number of threads.
        species_primer_candidates=candidates_dict,
        exclude_dict=round2_exclude_dict,
    )

# See PyCharm help at https://www.jetbrains.com/help/pycharm/
