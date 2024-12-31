
import os
from olivar.pipeline_review_redesign_getPrimer import build, tiling, panel_redesign_optimize, panel_save#, save, validate
import pickle
import re
import pandas as pd
import ast
import time
from datetime import datetime

from Bio import SeqIO
# import pandas as pd
# import os


species_primer_candidates = {}

# second round 避免第一轮的amp_1, 33, 56 区域
#round2_exclude_dict = {}
# round2_exclude_dict['amp_1'] = {'fp':(6592, 6615),'rp':()}
# round2_exclude_dict['amp_33'] = {'fp':(),'rp':(2155763,2155783)}
# round2_exclude_dict['amp_56'] = {'fp':(4326971,4326991),'rp':()}
# second round


# species_primer_candidates
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='引物设计代码')
    parser.add_argument('--first_design_file', help='the first_design_file')
    parser.add_argument('--bad_primer_path', help='bad_primer_dict path')
    parser.add_argument('--file_name', help='the file_name')
    args = parser.parse_args()
    first_design_file = args.first_design_file
    bad_primer_path = args.bad_primer_path
    file_name = args.file_name

    #config parameter
    title = file_name
    max_amp_len = 270  # Maximum amplicon length [420].
    min_amp_len = 210  # Minimum amplicon length. 0.9*{max-amp-len} if set as None. Minimum 120.
    w_egc = 1  # Weight for extreme GC content [1.0].
    w_lc = 1  # Weight for low sequence complexity [1.0].
    w_ns = 1  # Weight for non-specificity [1.0].
    w_var = 1  # Weight for variations [1.0]. 1,

    temperature = 60  # PCR annealing temperature [60.0].
    salinity = 0.18  # Concentration of monovalent ions in units of molar [0.18].
    dG_max = -10.8  # Maximum free energy change of a primer in kcal/mol [-11.8].  # dG_max=-11.8(min_len=15)
    min_GC = 0.4  # Minimum GC content of a primer [0.2].
    max_GC = 0.62  # Maximum GC content of a primer [0.75].
    min_complexity = 0.1  # Minimum sequence complexity of a primer [0.4].
    max_len = 29  # Maximum length of a primer [36].
    # check_var=True,  # Filter out primer candidates with variations within 5nt of 3' end [False].
    # Setting check_var=True is not recommended when a lot of variations are provided,
    # since this would significantly reduce the number of primer candidates.

    fP_prefix = ''  # Prefix of forward primer. Empty string '' by default.
    rP_prefix = ''  # Prefix of reverse primer. Empty string '' by default.

    seed = 10  # Random seed for optimizing primer design regions and primer dimer [10].
    threads = 1  # Number of threads. default 10
    #config parameter
    now = datetime.now()
    date_time = now.strftime("%Y-%m-%d, %H:%M:%S")
    output_panel = '/reference_data/primer_result/output/'
    if not os.path.exists(output_panel):
        os.mkdir(output_panel)
    # previous_design_path = './re_build/output/2024-10-29, 17:41:26/' # agent customized  #customized_panel_1.olvd
    # output_panel = os.path.join(previous_design_path,'redesign get primer 2')
    # if not os.path.exists(output_panel):
    #     os.mkdir(output_panel)
    all_plex_info_primer_dict = {}
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
        'min_complexity': min_complexity, 
        'max_len': max_len, 
        'fP_prefix': fP_prefix, 
        'rP_prefix': rP_prefix, 
        'seed': seed, 
        'threads': threads
    }
    pre_index = 0
    optimize_path = first_design_file
    # optimize_path = os.path.join(previous_design_path,'copy_2.olvd')
    with open(optimize_path,  'rb') as f:
        optimize_plex_info = pickle.load(f)


    redesign_amp_df = pd.read_csv(bad_primer_path)
    redesign_df = redesign_amp_df[:3]
    redesign_amp_dict = {}
    for k , v in zip(list(redesign_df['amp_id']),list(redesign_df['primer_type'])):
        redesign_amp_dict[k] = v
    # required parameters: plex names to be redesigned
    # redesign_plex_dict = {'amp_1': 'f', 'amp_2': 'fr'}
    redesign_plex_dict = {'amp_11': 'fr', 'amp_12': 'fr', 'amp_13': 'fr', 'amp_22': 'fr', 'amp_23': 'fr'}
    # redesign_plex_dict = redesign_amp_dict
    all_plex_info_optimize, learning_curve = panel_redesign_optimize(optimize_plex_info, redesign_plex_dict, config, species_primer_candidates, output_panel)
    df = panel_save(all_plex_info_optimize, config, species_primer_candidates,output_panel)

    # save design file
    design_out = {
        'config': config,
        'df': df, 
        'all_plex_info': all_plex_info_optimize, 
        'learning_curve': learning_curve
    }
    design_path = os.path.join(output_panel, f'{title}.pkl')
    with open(design_path, 'wb') as f:
        pickle.dump(design_out, f, protocol=4)
        print('Design file saved as %s' % design_path)
