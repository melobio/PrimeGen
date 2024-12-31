###### This is a sample Python script.

# Press ⌃R to execute it or replace it with your code.
# Press Double ⇧ to search everywhere for classes, files, tool windows, actions, and settings.
import os
os.chdir('./primer_agent')
from olivar.pipeline_review_redesign_dimer import build, tiling, panel_redesign_optimize, panel_save#, save, validate
import pickle
import re
import pandas as pd
import ast
import time
from datetime import datetime

from Bio import SeqIO
# import pandas as pd
# import os


# # get fasta list
# fasta_path = './hg38/out-Copy1.fasta'
# #fasta_path = './hg38/out_profound.fasta'

# fasta_list = {}
# for cds in SeqIO.parse(fasta_path, 'fasta'):
#     seq_raw = str(cds.seq)
    
#     fasta_list[cds.id] = seq_raw
    
        

# print('Get current directory',os.getcwd())

## species_primer_candidates
# tsv_path = './snp_target_typing/species'  #primer_agent/
# #tsv_path = './snp_mtb/20/primer3_files_20'  #primer_agent/
# tsv_names = os.listdir(tsv_path)

#print(os.path.abspath('./'))
species_primer_candidates = {}

# second round - avoid amp_1, 33, 56 regions from first round
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
#     df = pd.read_csv(t_path, sep="\t", header=None)  # ******* How to prevent first row from being column names *********
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

    #config parameter
    title = 'redesign_panel_dimer_1'
    max_amp_len=270  # Maximum amplicon length [420].
    min_amp_len=210 # Minimum amplicon length. 0.9*{max-amp-len} if set as None. Minimum 120.
    w_egc=1  # Weight for extreme GC content [1.0].
    w_lc=1 # Weight for low sequence complexity [1.0].
    w_ns=1 # Weight for non-specificity [1.0].
    w_var=1 # Weight for variations [1.0]. 1,

    temperature=60 # PCR annealing temperature [60.0].
    salinity=0.18 # Concentration of monovalent ions in units of molar [0.18].
    dG_max=-10.8  # Maximum free energy change of a primer in kcal/mol [-11.8].  # dG_max=-11.8(min_len=15) 
    min_GC=0.4 # Minimum GC content of a primer [0.2].
    max_GC=0.62 # Maximum GC content of a primer [0.75].
    min_complexity=0.1  # Minimum sequence complexity of a primer [0.4].
    max_len=29 # Maximum length of a primer [36].
    #check_var=True,  # Filter out primer candidates with variations within 5nt of 3' end [False].
    # Setting check_var=True is not recommended when a lot of variations are provided,
    # since this would significantly reduce the number of primer candidates.

    fP_prefix=''  # Prefix of forward primer. Empty string '' by default.
    rP_prefix=''  # Prefix of reverse primer. Empty string '' by default.

    seed=10  # Random seed for optimizing primer design regions and primer dimer [10].
    threads=1  # Number of threads. default 10
    #config parameter


    now = datetime.now()
    date_time = now.strftime("%Y-%m-%d, %H:%M:%S")
    #build fasta_list one by one
    # output_build = os.path.join('./re_build/build/',date_time)
    # os.mkdir(output_build)

    previous_design_path = './re_build/output/2024-10-29, 17:41:26/' # agent customized  #customized_panel_1.olvd
    output_panel = os.path.join(previous_design_path,'redesign')
    if not os.path.exists(output_panel):

        os.mkdir(output_panel)

    #
    all_plex_info_primer_dict = {}
    # #
    config = {
        #'ref_path': ref_path,  #
        #'out_path': out_path,  #
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
        'min_complexity': min_complexity, 
        'max_len': max_len, 
        'fP_prefix': fP_prefix, 
        'rP_prefix': rP_prefix, 
        'seed': seed, 
        'threads': threads
    }

    #
    pre_index = 0
    #
    

    #for i,(id,fa) in enumerate(fasta_list.items()):

        # build(
        #     #
        #     merge_id = id,

        #     fasta = fa, # Path to the fasta reference sequence.
        #     var_path = None,
        #     #var_path = './build_file/mtb_typing/20241011_var_plus1.csv', # Optional, path to the csv file of SNP coordinates and frequencies.
        #     # Required columns: "START", "STOP", "FREQ". "FREQ" is considered as 1.0 if empty. Coordinates are 1-based.
        #     BLAST_db = None,
        #     ## Must write the database name output from makeblastdb
        #     #BLAST_db = '/Users/wangyi/PycharmProjects/pcr_1/twt_设计/dt1/dt1_db/dt1_db', ## Must write the database name output from makeblastdb
        #     #Optional, path to the BLAST database.
        #     # Note that this path should end with the name of the BLAST database (e.g., "example_input/Human/GRCh38_primary").
        #     out_path = output_build,
        #     #out_path = './snp_target_typing', # Output directory.
        #     title = f'omim_{i}', # Name of the Olivar reference file [olivar-ref].
        #     threads = 1 # Number of threads.
        # )

        # print("debug")
        # #exclude_dict
        # exclude_dict = {}

        # #build_ref = os.path.join(output_build,f'omim_{i}.olvr')
        # #all_plex_info_primer,  cur_index = tiling(build_ref, exclude_dict, pre_index, config)

        # pre_index = cur_index
        #     #ref_path='./hsp65/hsp65-ref.olvr',  # Path to the Olivar reference file (.olvr).
        #     ref_path=os.path.join(output_build,f'omim_{i}.olvr'),
        #     #ref_path='./snp_mtb/20240110_mtb_snp_ref.olvr',
        #     #ref_path='./twt_设计/dt1/output/dt1-ref.olvr',

        #     #out_path='./hsp65/design',  # Output directory.
        #     #out_path='./snp_target_typing/output/20240226',
        #     #out_path='./twt_设计/dt1/output/',
        #     #out_path='./snp_target_typing/primer_3_get_primers',
        #     # os.mkdir(f'')
        #     # out_path=os.path.join('./re_build/output/',),
            

        #     # title='dt1-p-design',  # Name of this design [olivar-design].
        #     #title='snp_target_typing_species13_1tube_tm60_1',
        #     title='rebuild_p3_2',
        #     #title='dt1_ada-p-design',

        #     max_amp_len=270,  # Maximum amplicon length [420].
        #     min_amp_len=210,  # Minimum amplicon length. 0.9*{max-amp-len} if set as None. Minimum 120.
        #     w_egc=1,  # Weight for extreme GC content [1.0].
        #     w_lc=1,  # Weight for low sequence complexity [1.0].
        #     w_ns=1,  # Weight for non-specificity [1.0].
        #     w_var=1,  # Weight for variations [1.0]. 1,

        #     temperature=53,  # PCR annealing temperature [60.0].
        #     salinity=0.18,  # Concentration of monovalent ions in units of molar [0.18].
        #     dG_max=-11.8,  # Maximum free energy change of a primer in kcal/mol [-11.8]. 
        #     min_GC=0.4,  # Minimum GC content of a primer [0.2].
        #     max_GC=0.65,  # Maximum GC content of a primer [0.75].
        #     min_complexity=0.1,  # Minimum sequence complexity of a primer [0.4].
        #     max_len=23,  # Maximum length of a primer [36].
        #     #check_var=True,  # Filter out primer candidates with variations within 5nt of 3' end [False].
        #     # Setting check_var=True is not recommended when a lot of variations are provided,
        #     # since this would significantly reduce the number of primer candidates.

        #     fP_prefix='',  # Prefix of forward primer. Empty string '' by default.
        #     rP_prefix='',  # Prefix of reverse primer. Empty string '' by default.

        #     seed=10,  # Random seed for optimizing primer design regions and primer dimer [10].
        #     threads=10,  # Number of threads.
        # #     # species_primer_candidates=candidates_dict,
        # #     # exclude_dict=round2_exclude_dict,
        # # 

        #
    #all_plex_info_primer_dict.update(all_plex_info_primer)   #(all_plex_info_primer)
        # #
    
    #save candidates for redesign
    # candidates_path = os.path.join(output_panel, '%s.PGenC' % config['title'])
    # with open(candidates_path, 'wb') as f:
    #     pickle.dump(all_plex_info_primer_dict, f, protocol=4)
    #     print('Candidates file saved as %s' % candidates_path)
    
    #optimize(all_plex_info, config, species_candidates):
    #panel_optimize(all_plex_info, config)
    
    optimize_path = os.path.join(previous_design_path,'copy_2.olvd')
    with open(optimize_path,  'rb') as f:
        optimize_plex_info = pickle.load(f)
        
    # required parameters: plex names to be redesigned
    redesign_plex_dict = {'amp_1':'f','amp_2':'fr'}

    all_plex_info_optimize, learning_curve = panel_redesign_optimize(optimize_plex_info, redesign_plex_dict, config, species_primer_candidates, output_panel)
    df = panel_save(all_plex_info_optimize, config, species_primer_candidates,output_panel)

    # save design file
    design_out = {
        'config': config, 
        #'risk_arr': risk_arr, 
        #'gc_arr': gc_arr, 
        #'comp_arr': comp_arr, 
        #'hits_arr': hits_arr, 
        #'var_arr': var_arr, 
        #'all_loss': all_loss, 
        #'seq': seq_raw, 
        'df': df, 
        'all_plex_info': all_plex_info_optimize, 
        'learning_curve': learning_curve
    }
    design_path = os.path.join(output_panel, '%s.PGenRedesign' % config['title'])
    with open(design_path, 'wb') as f:
        pickle.dump(design_out, f, protocol=4)
        print('Design file saved as %s' % design_path)

# See PyCharm help at https://www.jetbrains.com/help/pycharm/
