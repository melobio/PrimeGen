###### This is a sample Python script.

# Press ⌃R to execute it or replace it with your code.
# Press Double ⇧ to search everywhere for classes, files, tool windows, actions, and settings.
import os
#os.chdir('./primer_agent')
from olivar.pipeline_review_snp import build, tiling#, save, validate

import re
import pandas as pd
import ast
import time
from datetime import datetime

## species_primer_candidates
# tsv_path = './snp_target_typing/species'  #primer_agent/
# #tsv_path = './snp_mtb/20/primer3_files_20'  #primer_agent/
# tsv_names = os.listdir(tsv_path)

#print(os.path.abspath('./'))
candidates_dict = {}

# second round 
round2_exclude_dict = {}
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

    now = datetime.now()
    date_time = now.strftime("%Y-%m-%d, %H:%M:%S")
    #build fasta_list one by one
    output_build = os.path.join('./re_build/cov19_build/',date_time)
    os.mkdir(output_build)
    output_panel = os.path.join('./re_build/cov19_output/',date_time)
    os.mkdir(output_panel)

    #build
    build(
        fasta_path = './cov19/NC_045512.2.fasta', # Path to the fasta reference sequence.
        #var_path = None,
        var_path = './cov19/covid19_del_1024_freq01.csv', # Optional, path to the csv file of SNP coordinates and frequencies.
        # Required columns: "START", "STOP", "FREQ". "FREQ" is considered as 1.0 if empty. Coordinates are 1-based.
        BLAST_db = None,
       
        #BLAST_db = '/Users/wangyi/PycharmProjects/pcr_1/twt_设计/dt1/dt1_db/dt1_db', ## 一定要写makeblastdb得出的数据库name
        #Optional, path to the BLAST database.
        # Note that this path should end with the name of the BLAST database (e.g., "example_input/Human/GRCh38_primary").
        out_path = output_build,
        #out_path = './snp_target_typing', # Output directory.
        title = 'cov19_snp', # Name of the Olivar reference file [olivar-ref].
        threads = 1 # Number of threads.
    )


    tiling(
        #ref_path='./hsp65/hsp65-ref.olvr',  # Path to the Olivar reference file (.olvr).
        ref_path = os.path.join(output_build,f'cov19_snp.olvr'),
        #ref_path='./snp_mtb/20240110_mtb_snp_ref.olvr',
        #ref_path='./twt_设计/dt1/output/dt1-ref.olvr',

        #out_path='./hsp65/design',  # Output directory.
        #out_path='./snp_target_typing/output/20240226',
        #out_path='./twt_设计/dt1/output/',
        #out_path='./snp_target_typing/primer_3_get_primers',
        out_path = output_panel,

        # title='dt1-p-design',  # Name of this design [olivar-design].
        #title='snp_target_typing_species13_1tube_tm60_1',
        title='snp_primer3',
        #title='dt1_ada-p-design',

        max_amp_len=370,  # Maximum amplicon length [420].
        min_amp_len=330,  # Minimum amplicon length. 0.9*{max-amp-len} if set as None. Minimum 120.
        w_egc=1,  # Weight for extreme GC content [1.0].
        w_lc=1,  # Weight for low sequence complexity [1.0].
        w_ns=1,  # Weight for non-specificity [1.0].
        w_var=100,  # Weight for variations [1.0]. 1,

        temperature=60,  # PCR annealing temperature [60.0].
        salinity=0.18,  # Concentration of monovalent ions in units of molar [0.18].
        dG_max=-10.8,  # Maximum free energy change of a primer in kcal/mol [-11.8]. 
        min_GC=0.4,  # Minimum GC content of a primer [0.2].
        max_GC=0.62,  # Maximum GC content of a primer [0.75].
        min_complexity=0.1,  # Minimum sequence complexity of a primer [0.4].
        max_len=23,  # Maximum length of a primer [36].
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