

from olivar import build, tiling#, save, validate
import os
import re
import pandas as pd
import ast

candidates_dict = {}

# second round
round2_exclude_dict = {}



if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='引物设计代码')
    parser.add_argument('--fasta_path', help='the path of fasta file')
    parser.add_argument('--out_path', help='the path of output')
    parser.add_argument('--gene_name', help='the gene name')
    parser.add_argument('--start_pos', type=int, nargs='?', const=1, default=0,help='the start_pos')
    parser.add_argument('--end_pos', type=int, nargs='?', const=1, default=0, help='the end_pos')
    parser.add_argument('--snp_path', help='the path of snp file',default=None)
    parser.add_argument('--max_amp_len', help='Maximum amplicon length',default=300)
    parser.add_argument('--min_amp_len', help='Minimum amplicon length',default=250)
    parser.add_argument('--temperature', help='PCR annealing temperature',default=65)
    parser.add_argument('--min_GC', help='Minimum GC content of a primer',default=0.4)
    parser.add_argument('--max_GC', help='Maximum GC content of a primer',default=0.6)
    parser.add_argument('--max_len', help='Maximum length of a primer',default=23)
    args = parser.parse_args()
    fasta_path =args.fasta_path
    out_path = args.out_path
    var_path = args.snp_path
    gene_name = args.gene_name
    if args.start_pos:
        build(
            fasta_path = fasta_path,
            var_path = var_path,
            BLAST_db = None,
            out_path = out_path,
            title = f'primer_design_{gene_name}', # Name of the Olivar reference file [olivar-ref].
            threads = 1 ,# Number of threads.
            start_pos=int(args.start_pos),
            end_pos=int(args.end_pos),
        )
    else:
        build(
            fasta_path=fasta_path,
            var_path=var_path,
            BLAST_db=None,
            out_path=out_path,
            title=f'primer_design_{gene_name}',  # Name of the Olivar reference file [olivar-ref].
            threads=1,  # Number of threads.
        )

    tiling(

        ref_path=f'/reference_data/primer_result/primer_design_{gene_name}.olvr',
        out_path=out_path,
        #title='mtb23_no_species-p-design',
        title=f'primer_design_{gene_name}_result',
        max_amp_len=int(args.max_amp_len),  # Maximum amplicon length [420].
        min_amp_len=int(args.min_amp_len),  # Minimum amplicon length. 0.9*{max-amp-len} if set as None. Minimum 120.
        w_egc=1,  # Weight for extreme GC content [1.0].
        w_lc=1,  # Weight for low sequence complexity [1.0].
        w_ns=2,  # Weight for non-specificity [1.0]. dt7k:5;
        w_var=1,  # Weight for variations [1.0].
        temperature=int(args.temperature),  # PCR annealing temperature [60.0].
        salinity=0.18,  # Concentration of monovalent ions in units of molar [0.18].
        dG_max=-11.8,  # Maximum free energy change of a primer in kcal/mol [-11.8].
        min_GC=int(args.min_GC)/100,  # Minimum GC content of a primer [0.2].
        max_GC=int(args.max_GC)/100,  # Maximum GC content of a primer [0.75].
        min_complexity=0.6,  # Minimum sequence complexity of a primer [0.4].
        max_len=int(args.max_len),  # Maximum length of a primer [36].
        fP_prefix='',  # Prefix of forward primer. Empty string '' by default.
        rP_prefix='',  # Prefix of reverse primer. Empty string '' by default.
        seed=10,  # Random seed for optimizing primer design regions and primer dimer [10].
        threads=10,  # Number of threads.
        species_primer_candidates=candidates_dict,
        exclude_dict=round2_exclude_dict,
    )

