import os
#os.chdir('./primer_agent')
from olivar.pipeline_review_1 import build, tiling, panel_optimize, panel_save#, save, validate
import pickle
import re
import pandas as pd
import ast
import time
from datetime import datetime

from Bio import SeqIO

#print(os.path.abspath('./'))
species_primer_candidates = {}


#round2_exclude_dict = {}
# round2_exclude_dict['amp_1'] = {'fp':(6592, 6615),'rp':()}
# round2_exclude_dict['amp_33'] = {'fp':(),'rp':(2155763,2155783)}
# round2_exclude_dict['amp_56'] = {'fp':(4326971,4326991),'rp':()}
# second round

# for t in tsv_names:
#     if t.startswith('.'):
#         continue
#     list_f_r = []
#     #t_path = os.path.join(os.path.abspath('./'),'snp_mtb','20','primer3_files_20',t) #'primer_agent',
#     t_path = os.path.join(os.path.abspath('./'),'snp_target_typing','species',t)
#     amp = t.split('_cds')[0]
#     print(amp)
#     df = pd.read_csv(t_path, sep="\t", header=None)  # *******如何让第一行数据不是列名字*********
#     df.columns = ['name', 'content']
#     primer_f_string = df[df['name'] == 'PRIMER_LEFT']['content'].values[0]
#     res_f = ast.literal_eval(primer_f_string)

#     # for f in res_f:
#     #     list_f.append(f['SEQUENCE'])
#     primer_r_string = df[df['name'] == 'PRIMER_RIGHT']['content'].values[0]
#     res_r = ast.literal_eval(primer_r_string)
#     # for r in res_r:
#     #     list_r.append(r['SEQUENCE'])
#     assert len(res_f)==len(res_r)
#     for i in range(len(res_f)):
#         list_f_r.append({'f':res_f[i]['SEQUENCE'].lower(), 'r':res_r[i]['SEQUENCE'].lower()})
#     candidates_dict[amp] = list_f_r
# print(len(candidates_dict))
# #print(candidates_dict)

# species_primer_candidates
if __name__ == "__main__":

    import argparse
    parser = argparse.ArgumentParser(description='引物设计代码')
    parser.add_argument('--fasta_path', help='the path of fasta file')
    parser.add_argument('--file_name', help='the file_name')
    args = parser.parse_args()
    fasta_path = args.fasta_path
    file_name = args.file_name
    # out_path = args.out_path

    # fasta_path = './hg38/out-Copy2.fasta'
    # fasta_path = './hg38/out_profound.fasta'

    fasta_list = {}
    for cds in SeqIO.parse(fasta_path, 'fasta'):
        seq_raw = str(cds.seq)

        fasta_list[cds.id] = seq_raw
    #config parameter
    title = file_name
    max_amp_len=300  # Maximum amplicon length [420].
    min_amp_len=210 # Minimum amplicon length. 0.9*{max-amp-len} if set as None. Minimum 120.
    w_egc=1  # Weight for extreme GC content [1.0].
    w_lc=1 # Weight for low sequence complexity [1.0].
    w_ns=1 # Weight for non-specificity [1.0].
    w_var=0 # Weight for variations [1.0]. 1,

    temperature=60 # PCR annealing temperature [60.0].
    salinity=0.18 # Concentration of monovalent ions in units of molar [0.18].
    dG_max=-11.8  # Maximum free energy change of a primer in kcal/mol [-11.8].  # dG_max=-11.8(min_len=15) #G_max=-10.8(min_len=18)
    min_GC=0.4 # Minimum GC content of a primer [0.2].
    max_GC=0.65 # Maximum GC content of a primer [0.75].
    # min_complexity=0.1  # Minimum sequence complexity of a primer [0.4].
    max_len=36 # Maximum length of a primer [36].
    #check_var=True,  # Filter out primer candidates with variations within 5nt of 3' end [False].
    # Setting check_var=True is not recommended when a lot of variations are provided,
    # since this would significantly reduce the number of primer candidates.
    min_atcg = 0.1
    max_atcg = 0.4
    fP_prefix=''  # Prefix of forward primer. Empty string '' by default.
    rP_prefix=''  # Prefix of reverse primer. Empty string '' by default.
    seed=10  # Random seed for optimizing primer design regions and primer dimer [10].
    threads=30  # Number of threads. default 10
    now = datetime.now()
    date_time = now.strftime("%Y-%m-%d, %H:%M:%S")
    #build fasta_list one by one
    # output_build = os.path.join('/reference_data/primer_result/build/',date_time)
    # os.mkdir(output_build)
    # output_panel = os.path.join('/reference_data/primer_result/output/',date_time)
    # os.mkdir(output_panel)
    output_build = '/reference_data/primer_result/build/'
    if not os.path.exists(output_build):
        os.mkdir(output_build)
    output_panel = '/reference_data/primer_result/output/'
    if not os.path.exists(output_panel):
        os.mkdir(output_panel)

    all_plex_info_primer_dict = {}
    # #
    config = {
        # 'ref_path': ref_path,  #
        # 'out_path': out_path,  #
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
        # 'min_complexity': min_complexity,
        'max_len': max_len,
        'fP_prefix': fP_prefix,
        'rP_prefix': rP_prefix,
        'seed': seed,
        'threads': threads,
        'BLAST_db': '/reference_data/db/hg38'
    }
    pre_index = 0
    #

    for i, (id, fa) in enumerate(fasta_list.items()):
        build(
            #
            config,
            merge_id=id,

            fasta=fa,  # Path to the fasta reference sequence.
            var_path=None,
            # var_path = './build_file/mtb_typing/20241011_var_plus1.csv', # Optional, path to the csv file of SNP coordinates and frequencies.
            # Required columns: "START", "STOP", "FREQ". "FREQ" is considered as 1.0 if empty. Coordinates are 1-based.
            BLAST_db=config['BLAST_db'],

            # BLAST_db = '/Users/wangyi/PycharmProjects/pcr_1/twt_设计/dt1/dt1_db/dt1_db', ## 
            # Optional, path to the BLAST database.
            # Note that this path should end with the name of the BLAST database (e.g., "example_input/Human/GRCh38_primary").
            out_path=output_build,
            # out_path = './snp_target_typing', # Output directory.
            title=f'omim_{i}',  # Name of the Olivar reference file [olivar-ref].
            threads=config['threads']  # Number of threads.

        )

        print("debug")
        # exclude_dict
        exclude_dict = {}

        build_ref = os.path.join(output_build, f'omim_{i}.olvr')
        all_plex_info_primer, cur_index = tiling(build_ref, exclude_dict, pre_index, config)

        pre_index = cur_index

        all_plex_info_primer_dict.update(all_plex_info_primer)  # (all_plex_info_primer)

    all_plex_info_optimize, learning_curve = panel_optimize(all_plex_info_primer_dict, config,
                                                            output_panel)  # species_primer_candidates,
    df = panel_save(all_plex_info_optimize, config, output_panel)  # species_primer_candidates,

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

