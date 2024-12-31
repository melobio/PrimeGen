###### This is a sample Python script.

# Press ⌃R to execute it or replace it with your code.
# Press Double ⇧ to search everywhere for classes, files, tool windows, actions, and settings.
import os
# os.chdir('./primer_agent')
from olivar.pipeline_review_2 import build, tiling, panel_optimize, panel_save#, save, validate
import pickle
import re
import pandas as pd
import ast
import time
from datetime import datetime

from Bio import SeqIO
# import pandas as pd
# import os


if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser(description='引物设计代码')
    parser.add_argument('--fasta_path', help='the path of fasta file')
    parser.add_argument('--snp_file_path', help='the path of snp file')
    parser.add_argument('--file_name', help='the file_name')
    parser.add_argument('--output_path', help='the output_path')
    args = parser.parse_args()
    fasta_path = args.fasta_path
    file_name = args.file_name
    output_path = args.output_path
    snp_file_path = args.snp_file_path
    fasta_list = {}
    for cds in SeqIO.parse(fasta_path, 'fasta'):
        seq_raw = str(cds.seq)

        fasta_list[cds.id] = seq_raw

    species_primer_candidates = {}
    #config parameter
    title = file_name
    max_amp_len=370  # Maximum amplicon length [420].
    min_amp_len=330 # Minimum amplicon length. 0.9*{max-amp-len} if set as None. Minimum 120.
    w_egc=1  # Weight for extreme GC content [1.0].
    w_lc=1 # Weight for low sequence complexity [1.0].
    w_ns=1 # Weight for non-specificity [1.0].
    w_var=1 # Weight for variations [1.0]. 1,

    temperature=60 # PCR annealing temperature [60.0].
    salinity=0.18 # Concentration of monovalent ions in units of molar [0.18].
    dG_max=-11.8  # Maximum free energy change of a primer in kcal/mol [-11.8].  # dG_max=-11.8(min_len=15) #G_max=-10.8(min_len=18)
    min_GC=0.4 # Minimum GC content of a primer [0.2].
    max_GC=0.62 # Maximum GC content of a primer [0.75].

    #
    min_atcg = 0.1
    max_atcg = 0.4
    #

    #min_complexity=0.4  # Minimum sequence complexity of a primer [0.4].
    max_len=36 # Maximum length of a primer [36].
    #check_var=True,  # Filter out primer candidates with variations within 5nt of 3' end [False].
    # Setting check_var=True is not recommended when a lot of variations are provided,
    # since this would significantly reduce the number of primer candidates.

    fP_prefix=''  # Prefix of forward primer. Empty string '' by default.
    rP_prefix=''  # Prefix of reverse primer. Empty string '' by default.

    seed=10  # Random seed for optimizing primer design regions and primer dimer [10].
    threads=1  # Number of threads. default 10
    #config parameter
    # use_snp_cluster
    use_snp_cluster = False
    # use_snp_cluster


    now = datetime.now()
    output_build = '/reference_data/primer_result/build/'
    if not os.path.exists(output_build):
        os.mkdir(output_build)
    output_panel = output_path
    if not os.path.exists(output_panel):
        os.mkdir(output_panel)
    # date_time = now.strftime("%Y-%m-%d, %H:%M:%S")
    #build fasta_list one by one
    # output_build = os.path.join('./re_build/cov19_build/',date_time)
    # os.mkdir(output_build)
    # output_panel = os.path.join('./re_build/cov19_output/whole_panel_redesign/',date_time)
    # os.mkdir(output_panel)
    all_plex_info_primer_dict = {}
    # #
    config = {
        'title': title,  
        'max_amp_len': max_amp_len, 
        'min_amp_len': min_amp_len, 
        'w_egc': w_egc, 
        'w_lc': w_lc, 
        'w_ns': w_ns, 
        'w_var': w_var, 
        'temperature': temperature, 
        'salinity': salinity, 
        'dG_max': dG_max, 
        'min_GC': min_GC, 
        'max_GC': max_GC, 
        'min_atcg': min_atcg,
        'max_atcg': max_atcg,
        #'min_complexity': min_complexity, 
        'max_len': max_len, 
        'fP_prefix': fP_prefix, 
        'rP_prefix': rP_prefix, 
        'seed': seed, 
        'threads': threads,
        'use_snp_cluster': use_snp_cluster
    }

    pre_index = 0
    for i,(id,fa) in enumerate(fasta_list.items()):

        build(
            #
            config,
            merge_id = id,

            fasta = fa, # Path to the fasta reference sequence.
            var_path = snp_file_path,
            #var_path = './build_file/mtb_typing/20241011_var_plus1.csv', # Optional, path to the csv file of SNP coordinates and frequencies.
            # Required columns: "START", "STOP", "FREQ". "FREQ" is considered as 1.0 if empty. Coordinates are 1-based.
            BLAST_db = None,
         
            #BLAST_db = '/Users/wangyi/PycharmProjects/pcr_1/twt_设计/dt1/dt1_db/dt1_db', 
            #Optional, path to the BLAST database.
            # Note that this path should end with the name of the BLAST database (e.g., "example_input/Human/GRCh38_primary").
            out_path = output_build,
            #out_path = './snp_target_typing', # Output directory.
            title = f'cov_whole_panel_{i}', # Name of the Olivar reference file [olivar-ref].
            threads = 1 # Number of threads.
        )
        print("debug")
        #exclude_dict
        exclude_dict = {}
        build_ref = os.path.join(output_build,f'cov_whole_panel_{i}.olvr')
        all_plex_info_primer,  cur_index = tiling(build_ref, exclude_dict, pre_index, config)
        pre_index = cur_index
        all_plex_info_primer_dict.update(all_plex_info_primer)   #(all_plex_info_primer)
        # #
    all_plex_info_optimize, learning_curve = panel_optimize(all_plex_info_primer_dict, config, output_panel) #species_primer_candidates,
    df = panel_save(all_plex_info_optimize, config, output_panel) #  species_primer_candidates,

    # save design file
    design_out = {
        'config': config,
        'df': df, 
        'all_plex_info': all_plex_info_optimize, 
        'learning_curve': learning_curve
    }
    design_path = os.path.join(output_panel, '%s.olvd' % config['title'])
    with open(design_path, 'wb') as f:
        pickle.dump(design_out, f, protocol=4)
        print('Design file saved as %s' % design_path)

