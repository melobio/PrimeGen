#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Main workflow of Olivar tiling.
Architecture:
main.py
    build()
    tiling()
        design_context_seq()
            generate_context()
                find_min_loc()
        get_primer()
        optimize()
        save()
    visualize()
    validate()
'''

import os
CURRENT_DIR = os.path.dirname(__file__)

import json
import pickle
import random
import copy
import multiprocessing
from time import time
from math import ceil, floor
from copy import deepcopy

import pandas as pd
import numpy as np
from numpy.random import rand, default_rng
from Bio import SeqIO
from matplotlib import pyplot as plt
from matplotlib.widgets import Slider

from . import basic
from .  import design#import design
from .ncbi_tools import BLAST_batch_short, ns_simulation
from .design import PrimerSetBadnessFast

#plt.rcParams.update({'font.family':'Arial'})

PRIMER_DESIGN_LEN = 40 # length of primer design region [40]
UPPER_GC = 0.75 # extreme GC upper bond [0.75]
LOWER_GC = 0.25 # extreme GC lower bond [0.75]
LOW_COMPLEXITY = 0.4 # low complexity lower bond [0.4]
CHOICE_RATE = 0.3 # bottom CHOICE_RATE lowest risk primer design regions are choosen from candidate region
RISK_TH = 0.1 # top RISK_TH highest risk primer design regions are considered as loss

CLUSTER_GAP = 300  # snp clustering gap

## ada
def index_to_seq(d_l):
    #d_l = list(index_dict.values())
    d_s = ["{:03}".format(i) for i in d_l]
    res = "".join(d_s)
    return res

def seq_to_index(seq):
    res_l = []
    if len(seq)%3!=0:
        raise ValueError('seq length should be divisible by 3.')
    else:
        for i in range(0, len(seq), 3):
            res_l.append(int(seq[i:i+3]))
            #print(seq[i:i+2])
    return res_l

def generate_new_combination(d,mu,dd):
    m = np.array(list(dd.values())).mean()
    mutant = []
    for k,v in d.items():
        if random.random() < mu * (dd[k]/m):
            mutant.append(floor(rand() * dd[k]))
        else:
            mutant.append(v)
    return mutant

def generate_new_combination_2(d,mu,dd):
    #m = np.array(list(dd.values())).mean()
    mutant = []
    c_tube = list(d.keys())
    mutplex = random.choice(c_tube)
    x = c_tube.index(mutplex)
    mutindex = floor(rand() * dd[mutplex])
    res = list(d.values())
    res[x] = mutindex
    return res

def compute_badness(all_plex_info, curr_tube, r, existing_primer):
    _curr_fp = {}
    _curr_rp = {}
    _curr_badness = 0
    root_ind = seq_to_index(r)
    # roots_inds.append(root_ind)
    if len(root_ind) != len(curr_tube):
        raise ValueError('plex index 出错.')
    # new_index = dict(zip(list(curr_index.keys()),root_ind))  # # #
    for plex_id, ind in zip(curr_tube, root_ind):
        _fp = all_plex_info[plex_id]['optimize'][ind][0]
        _rp = all_plex_info[plex_id]['optimize'][ind][1]
        _curr_fp[plex_id] = _fp['seq']
        _curr_rp[plex_id] = _rp['seq']
        _curr_badness += _fp['badness'] + _rp['badness']

        #cost += 1

    _inter_badness, _comp_badness = PrimerSetBadnessFast(list(_curr_fp.values()),
                                                         list(_curr_rp.values()), existing_primer)
    _curr_badness += _inter_badness
    return _curr_badness, _comp_badness


def _recombine_population(gen, recomb_rate):
    # If only one member of population, can't do any recombining
    if len(gen) == 1:
        return gen

    random.shuffle(gen)
    ret = []
    for i in range(0, len(gen) - 1, 2):
        strA = []
        strB = []
        switch = False

        for ind in range(len(gen[i])):
            if random.random() < recomb_rate:
                switch = not switch

            # putting together recombinants
            if switch:
                strA.append(gen[i][ind])
                strB.append(gen[i + 1][ind])
            else:
                strB.append(gen[i][ind])
                strA.append(gen[i + 1][ind])

        ret.append("".join(strA))
        ret.append("".join(strB))
    return ret
def _recombine_population_c(gen, recomb_rate):
    # If only one member of population, can't do any recombining
    if len(gen) == 1:
        return gen

    random.shuffle(gen)
    ret = []
    for i in range(0, len(gen) - 1, 2):
        strA = []
        strB = []
        switch = False
        #
        gen[i] = seq_to_index(gen[i])
        gen[i + 1] = seq_to_index(gen[i + 1])
        #
        
        for ind in range(len(gen[i])):
            if random.random() < recomb_rate:
                switch = not switch

            # putting together recombinants
            if switch:
                strA.append(gen[i][ind])
                strB.append(gen[i + 1][ind])
            else:
                strB.append(gen[i][ind])
                strA.append(gen[i + 1][ind])

        #
        ret.append("".join(index_to_seq(strA)))
        ret.append("".join(index_to_seq(strB)))
    return ret
## ada

def build(fasta_path, var_path, BLAST_db, out_path, title, threads):
    '''
    Build the Olivar reference file for tiled amplicon design
    Input:
        fasta_path: Path to the fasta reference sequence.
        var_path: Optional, path to the csv file of SNP coordinates and frequencies. 
            Required columns: "START", "STOP", "FREQ". "FREQ" is considered as 1.0 if empty. Coordinates are 1-based.
        BLAST_db: Optional, path to the BLAST database. 
            Note that this path should end with the name of the BLAST database (e.g., "example_input/Human/GRCh38_primary").
        out_path: Output directory [./]. 
        title: Name of the Olivar reference file [olivar-ref]. 
        threads: Number of threads. 
    '''
    n_cpu = threads

    if not os.path.exists(out_path):
        os.makedirs(out_path)

    # all coordinates are 1-based and inclusive
    word_size = 28
    offset = 14 # word_size should be divisible by offset
    print('building Olivar reference with word_size=%d, offset=%d' % (word_size, offset))

    n_cycle = word_size/offset
    if n_cycle != int(n_cycle):
        raise ValueError('word_size should be divisible by offset.')
    else:
        n_cycle = int(n_cycle)

    for record in SeqIO.parse(fasta_path, 'fasta'):
        seq_raw = str(record.seq)

    # load SNP/iSNV coordinates
    var_arr = np.zeros(len(seq_raw))
    seq_raw = seq_raw.lower() # SNPs in upper case
    if var_path:
        seq_raw = list(seq_raw)
        df = pd.read_csv(var_path, sep=',')
        for i, row in df.iterrows():
            if np.isnan(row['FREQ']):
                freq = 1.0
            else:
                freq = row['FREQ']**0.5
            var_arr[int(row['START'])-1:int(row['STOP'])] += freq
            # capitalize SNPs
            for pos in range(int(row['START'])-1, int(row['STOP'])):
                seq_raw[pos] = seq_raw[pos].upper()
        seq_raw = ''.join(seq_raw)
    else:
        print('No sequence variation coordinates provided, skipped.')

    # calculate risk of GC, non-specificity, complexity
    seq_len_temp = (len(seq_raw)//offset) * offset # length of temporary sequence
    start_temp = (len(seq_raw) - seq_len_temp)//2 + 1 # start coordinate of sequence to be processed as words
    stop_temp = start_temp + seq_len_temp - 1

    # fetch all words
    all_word = []
    for pos in range(0, seq_len_temp-word_size+1, offset):
        word_start = start_temp+pos-1
        word = seq_raw[word_start : word_start+word_size]
        all_word.append(word)

    # get GC, complexity, BLAST hits of each word
    with multiprocessing.Pool(processes=n_cpu) as pool:
        print('calculating GC content and sequence complexity...')
        all_gc = np.array(pool.map(basic.get_GC, all_word))
        all_complexity = np.array(pool.map(basic.get_complexity, all_word))
    if BLAST_db:
        all_hits, _ = BLAST_batch_short(all_word, db=BLAST_db, n_cpu=n_cpu, mode='rough') # tabular output
    else:
        print('No BLAST database provided, skipped.')
        all_hits = np.zeros(len(all_word))

    # calculate score array
    start = start_temp + (n_cycle-1)*offset # start coordinate of sequence with score array
    stop = stop_temp - (n_cycle-1)*offset
    seq_len = stop - start + 1 # length of final sequence
    gc_arr = np.zeros(seq_len)
    comp_arr = np.zeros(seq_len)
    hits_arr = np.zeros(seq_len)
    for pos in range(seq_len):
        n = pos//offset
        gc_arr[pos] = np.sum(all_gc[n:n+n_cycle])
        comp_arr[pos] = np.sum(all_complexity[n:n+n_cycle])
        hits_arr[pos] = np.sum(all_hits[n:n+n_cycle])
    gc_arr = gc_arr/n_cycle
    hits_arr = hits_arr/n_cycle

    olv_ref = {
        'seq': seq_raw, 
        'start': start, 
        'stop': stop, 
        'gc_arr': gc_arr, 
        'comp_arr': comp_arr, 
        'var_arr': var_arr[start-1:stop], 
        'hits_arr': hits_arr, 
        'all_hits': all_hits
    }
    save_path = os.path.join(out_path, '%s.olvr' % title)
    with open(save_path, 'wb') as f:
        pickle.dump(olv_ref, f, protocol=4)
        print('Reference file saved as %s' % save_path)


def find_min_loc(risk_arr, start, stop, rng):
    '''
    find the next primer design region within start and stop.
    '''
    score = []
    loc = np.array(range(start, stop-PRIMER_DESIGN_LEN+2))
    for i in loc:
        score.append(sum(risk_arr[i-1: i-1+PRIMER_DESIGN_LEN]))
    if len(loc) == 1:
        idx_rand = 0
    elif CHOICE_RATE < 1:
        k = int(len(loc)*CHOICE_RATE)
        if k == 0:
            k = 1
        idx_k_lowest = np.argpartition(score, k)[:k]
        idx_rand = rng.choice(idx_k_lowest)
    else:
        idx_rand = rng.choice(range(len(loc)))
    primer_start = loc[idx_rand]
    return primer_start, primer_start+PRIMER_DESIGN_LEN-1, score[idx_rand]


def generate_context(SOI):
    '''
    find a valid set of PDRs that cover the whole sequence of interest (SOI).
    '''
    a_start, a_stop, risk_arr, max_amp_len, min_amp_len, var_ar, seed = SOI #start, stop,
    rng = default_rng(seed)

    print(var_ar)
    #check_zero = np.array(var_ar)
    if np.all(var_ar==0):
        clu = []
        clu_list = np.where(var_ar==0)[0][a_start:a_stop]
        clu.append(clu_list)
        print(clu)
    else:
        clu = snp_clu(var_ar)  # snp cluster


    all_context_seq = [] # coordinates of all context sequences
    all_risk = []
    clu_mun = 1
    # cellecting context seq according to snp clustering result
    for sub_clu in clu:
        start = sub_clu[0]
        stop = sub_clu[-1]
        # cellecting context seq according to snp clustering result
        # generate the first pair of primer design region
        fp_start, fp_stop, fp_risk = find_min_loc(risk_arr, start, start+3*PRIMER_DESIGN_LEN-1, rng)
        rp_stop, rp_start, rp_risk = find_min_loc(risk_arr, fp_start+min_amp_len-PRIMER_DESIGN_LEN, fp_start+max_amp_len-1, rng)
        all_context_seq.append((fp_start, fp_stop, rp_stop, rp_start))
        all_risk.append((fp_risk, rp_risk))

        # generate the second pair of primer design region
        left_lim = max(fp_stop+1, rp_start+2*PRIMER_DESIGN_LEN-max_amp_len+1) # leave space for the next fP design region
        if left_lim > stop:
            print(f'cluster {clu_mun} needs only 1 amplicon ')
            continue
        fp_start, fp_stop, fp_risk = find_min_loc(risk_arr, left_lim, rp_stop-1, rng)
        left_lim = max(fp_start+min_amp_len-PRIMER_DESIGN_LEN, rp_start+PRIMER_DESIGN_LEN+1)
        right_lim = fp_start+max_amp_len-1
        if right_lim > stop:
            #raise ValueError('sequence provided is too short.')
            print(f'cluster {clu_mun} needs only 1 amplicon ')
            continue
        rp_stop, rp_start, rp_risk = find_min_loc(risk_arr, left_lim, right_lim, rng)
        all_context_seq.append((fp_start, fp_stop, rp_stop, rp_start))
        all_risk.append((fp_risk, rp_risk))

        # design the rest context sequences
        cs_id = 2 # id of the next context sequence to be designed (0-based)
        while True:
            # fp
            # all_context_seq[cs_id-2][3]: rp_start of the last last context sequence
            left_lim = max(all_context_seq[cs_id-2][3]+1, rp_start+2*PRIMER_DESIGN_LEN-max_amp_len+1)
            fp_start, fp_stop, fp_risk = find_min_loc(risk_arr, left_lim, rp_stop-1, rng)
            # rp
            left_lim = max(fp_start+min_amp_len-PRIMER_DESIGN_LEN, rp_start+PRIMER_DESIGN_LEN+1)
            right_lim = fp_start+max_amp_len-1
            if right_lim > stop:
                break
            rp_stop, rp_start, rp_risk = find_min_loc(risk_arr, left_lim, right_lim, rng)

            all_context_seq.append((fp_start, fp_stop, rp_stop, rp_start))
            all_risk.append((fp_risk, rp_risk))
            cs_id += 1

        clu_mun += 1

    all_context_seq = np.array(all_context_seq)
    all_risk = np.array(all_risk)
    all_risk_flatten = all_risk.flatten()
    k = int(len(all_risk_flatten) * (1-RISK_TH))
    loss = sum(np.partition(all_risk_flatten, k)[k:]**2)
    #loss = sum(all_risk_flatten[all_risk_flatten>=RISK_TH])
    return all_context_seq, all_risk, loss

def snp_clu(a_r):
    '''
    snp聚类
    '''
    cluster = []
    var_index = np.where(a_r!=0)[0]
    last_pos = a_r[0]
    tmp_clu = []
    tmp_clu.append(var_index[0])
    for i in range(1,len(var_index)):
        if (var_index[i] - var_index[i-1]) < CLUSTER_GAP:

            tmp_clu.append(var_index[i])
        else:
            cluster.append(tmp_clu)
            tmp_clu = []
            tmp_clu.append(var_index[i])
    #print(len(cluster))
    return cluster

def design_context_seq(config):
    '''
    PDR optimization. 
    '''
    ref_path = config['ref_path']
    max_amp_len = config['max_amp_len']
    w_gc = config['w_egc']
    w_complexity = config['w_lc']
    w_hits = config['w_ns']
    w_var = config['w_var']
    seed = config['seed']
    n_cpu = config['threads']
    
    # set random number generator
    rng_parent = default_rng(seed)
    
    with open(ref_path,  'rb') as f:
        olv_ref = pickle.load(f)

    seq_raw = olv_ref['seq']
    start = olv_ref['start'] # start of design region, 1-based, closed
    stop = olv_ref['stop'] # stop of design region, 1-based, closed
    gc_arr = olv_ref['gc_arr']
    comp_arr = olv_ref['comp_arr']
    var_arr = olv_ref['var_arr']
    hits_arr = olv_ref['hits_arr']
    print('successfully loaded reference file %s' % ref_path)
    
    if not config['min_amp_len']:
        min_amp_len = int(max_amp_len*0.9)
    else:
        min_amp_len = config['min_amp_len']
    if min_amp_len < 3*PRIMER_DESIGN_LEN:
        raise ValueError('min_amp_len is too small')

    gc_arr = np.logical_or(gc_arr<LOWER_GC, gc_arr>UPPER_GC).astype(int)
    comp_arr = (comp_arr < LOW_COMPLEXITY).astype(int)
    if max(hits_arr) != 0:
        hits_arr = hits_arr/max(hits_arr) # normalize

    # make arrays the same length as seq_raw
    gc_arr = w_gc*0.5*np.concatenate((np.zeros(start-1), gc_arr, np.zeros(len(seq_raw)-stop)))
    comp_arr = w_complexity*0.5*np.concatenate((np.zeros(start-1), comp_arr, np.zeros(len(seq_raw)-stop)))
    hits_arr = w_hits*np.concatenate((np.zeros(start-1), hits_arr, np.zeros(len(seq_raw)-stop)))
    var_arr = w_var*10*np.concatenate((np.zeros(start-1), var_arr, np.zeros(len(seq_raw)-stop)))

    # construct risk array (first row is risk, second row is coordinate on seq_rawy)
    risk_arr = gc_arr + comp_arr + hits_arr + var_arr

    #N = 500*len(risk_arr)//max_amp_len # number of primer sets to generate
    N = 100
    rand_int = rng_parent.integers(2**32, size=N) # random seeds for each iteration

    # single thread

    tik = time()
    with multiprocessing.Pool(processes=n_cpu) as pool:
        batch = [(start, stop, risk_arr, max_amp_len, min_amp_len, var_arr, rand_int[i]) for i in range(N)] #start, stop,
        design = pool.map(generate_context, batch)
    print('finished in %.3fs' % (time()-tik))

    # find the best arangement of primer design regions
    best_design = sorted(design, key=lambda x:x[2])
    all_context_seq, all_risk, loss = best_design[0]
    print('Best total risk: %.3f' % loss)

    # plot
    all_loss = [d[2] for d in design]
    all_loss.sort(reverse=True)
    fig = plt.figure(figsize=(8, 5))
    plt.plot(all_loss)
    plt.title('Optimization of primer design regions', fontsize=14)
    plt.xlabel('iterations (sorted)', fontsize=14)
    plt.ylabel('Loss', fontsize=14)
    save_path = os.path.join(config['out_path'], '%s.png' % config['title'])
    plt.savefig(save_path, bbox_inches='tight')
    plt.close()
    print('PDR optimization figure saved as %s' % save_path)
    # prepare output
    all_plex_info = {}
    for i, (context_seq, risk) in enumerate(zip(all_context_seq, all_risk)):
        cs = seq_raw[context_seq[0]-1:context_seq[3]] # actual context sequence
        plex_info = {
            #'tube': i%2 + 1,
            'tube': 1,
            'context_seq_coords': (context_seq[0], context_seq[3]),
            'primers_coords': tuple(context_seq), 
            'insert_coords': (context_seq[1]+1, context_seq[2]-1), 
            'context_seq': cs, 
            'context_seq_revcomp': basic.revcomp(cs), 
            'fP_design': seq_raw[context_seq[0]-1:context_seq[1]], 
            'rP_design': basic.revcomp(seq_raw[context_seq[2]-1:context_seq[3]]), 
            'risk': tuple(risk)
        }
        all_plex_info['amp_%d' % (i+1)] = plex_info
    print('total amplicons: %d' % (i+1))
    #print('total amplicons: %d' % (len(all_risk) + 1))
    cover_start = all_context_seq[0][1]+1
    cover_stop = all_context_seq[-1][2]-1
    print('covered region: %d:%d' % (cover_start, cover_stop))
    print('coverage of reference sequence: %.3f%%' % (100*(cover_stop-cover_start+1)/len(seq_raw)))
    #return seq_raw, risk_arr, all_context_seq, all_risk, all_plex_info
    return all_plex_info, risk_arr, gc_arr, comp_arr, hits_arr, var_arr, all_loss, seq_raw


def get_primer(all_plex_info, config):
    '''
    Generate candidates of fP and rP for each plex
    Input:
        all_plex_info: output of design_context_seq(all_gene_info, config), or load from 'all_plex_info.pkl'
        config: design configurations, load from 'config.json'
    Output:
        all_plex_info and corresponding 'all_plex_info_primer.pkl': add primer candidates to input all_plex_info
    '''
    print('Generating primer candidates...')

    temperature_fP = config['temperature']
    temperature_rP = config['temperature']
    salinity = config['salinity']
    fP_adapter = config['fP_prefix']
    rP_adapter = config['rP_prefix']
    default_fP_setting = {
        'dG_max': config['dG_max'], 
        'min_GC': config['min_GC'], 
        'max_GC': config['max_GC'], 
        'min_complexity': config['min_complexity'], 
        'max_len': config['max_len'], 
        'check_SNP': True
    }
    default_rP_setting = deepcopy(default_fP_setting)
    check_BLAST = False

    # generate primer candidates
    fP_generator = design.primer_generator(temperature_fP, salinity)
    rP_generator = design.primer_generator(temperature_rP, salinity)
    n = 0
    for plex_id, plex_info in all_plex_info.items():
        # set primer prefix
        fP_prefix = fP_adapter
        rP_prefix = rP_adapter
        n += 1
        
        # generate fP
        fP_design = plex_info['fP_design'] # sequence to generate fP candidates
        fP_setting = deepcopy(default_fP_setting)
        fP, fail = fP_generator.get(fP_design, fP_prefix, check_BLAST=check_BLAST, **fP_setting)
        while len(fP) == 0:
            print('Fail to generate fP, update setting: %s' % plex_id)
            # detect failure mode with the most cases
            fail_mode = max(fail, key=fail.get)
            # update setting
            if fail_mode == 'SNP_fail' and fP_setting['check_SNP'] == True:
                fP_setting['check_SNP'] = False
            elif fail_mode == 'GC_fail':
                fP_setting['min_GC'] -= 0.05
                fP_setting['max_GC'] += 0.05
            elif fail_mode == 'complexity_fail':
                fP_setting['min_complexity'] -= 0.05
            # generate primers again
            fP, fail = fP_generator.get(fP_design, fP_prefix, check_BLAST=check_BLAST, **fP_setting)
        plex_info['fP_candidate'] = fP
        plex_info['fP_setting'] = fP_setting
        
        # generate rP
        rP_design = plex_info['rP_design'] # sequence to generate rP candidates
        rP_setting = deepcopy(default_rP_setting)
        rP, fail = rP_generator.get(rP_design, rP_prefix, check_BLAST=check_BLAST, **rP_setting)
        while len(rP) == 0:
            print('Fail to generate rP, update setting: %s' % plex_id)
            # detect failure mode with the most cases
            fail_mode = max(fail, key=fail.get)
            # update setting
            if fail_mode == 'SNP_fail' and rP_setting['check_SNP'] == True:
                rP_setting['check_SNP'] = False
            elif fail_mode == 'GC_fail':
                rP_setting['min_GC'] -= 0.05
                rP_setting['max_GC'] += 0.05
            elif fail_mode == 'complexity_fail':
                rP_setting['min_complexity'] -= 0.05
            # generate primers again
            rP, fail = rP_generator.get(rP_design, rP_prefix, check_BLAST=check_BLAST, **rP_setting)
        plex_info['rP_candidate'] = rP
        plex_info['rP_setting'] = rP_setting
        
    # with open('all_plex_info_primer.pkl', 'wb') as f:
    #     pickle.dump(all_plex_info, f)
    print('Done\n')
    return all_plex_info


def optimize(all_plex_info, config):
    '''
    Optimize fP and rP with the SADDLE algorithm.
    Input:
        all_plex_info: output of get_primer(all_plex_info, config), or load from 'all_plex_info_primer.pkl'
        config: design configurations, load from 'config.json'
    Output:
        all_plex_info and corresponding 'all_plex_info_optimize.pkl': add optimized fP and rP to each plex of input all_plex_info
    '''
    print('Optimizing primer dimers...')
    #set random seed
    # np.random.seed(config['seed'])
    # random.seed(config['seed'])

    #n_tube = 2
    #ada parameter
    cost = 0
    cost_budge = 2000
    measured_index_list = []
    batch_size = 2000
    eval_batch_size = 20
    c_threshold = 2  ### 这个参数需要摸索
    top_num = 50
    recomb_rate = 0
    #ada parameter
    
    
    n_tube = 1
    pool_size = max(len(all_plex_info)/n_tube - 100, 0)
    InitSATemp = 1000 + 10*pool_size
    NUMSTEPS = 10 + int(pool_size/10)
    ZEROSTEPS = NUMSTEPS
    TimePerStep = 1000
    # InitSATemp = config['InitSATemp']
    # NUMSTEPS = config['NumSteps']
    # ZEROSTEPS = config['ZeroSteps']
    # TimePerStep = config['TimePerStep']

    # generate primer pair candidates (fP-rP)
    n_pair = []
    for plex_id, plex_info in all_plex_info.items():
        optimize = []
        plex_amp_len = []
        for i, fp in plex_info['fP_candidate'].iterrows():
            for j, rp in plex_info['rP_candidate'].iterrows():
                insert_len = plex_info['insert_coords'][1] - plex_info['insert_coords'][0] + 1
                # check amplicon length
                amp_len = fp['primer_len'] + fp['dist'] + insert_len + rp['dist'] + rp['primer_len']
                plex_amp_len.append(amp_len)
                optimize.append([fp, rp])
        plex_info['optimize'] = optimize
        n_pair.append(len(optimize))
    print('total primer pairs %d' % sum(n_pair))
    print('average pairs per plex %.2f' % (sum(n_pair)/len(n_pair)))
    
    # optimize each tube
    all_lc = []
    tik = time()
    for i_tube in range(1, n_tube+1):
        print('\npool %d, simulated annealing...' % i_tube)
        # plex_id of current tube
        curr_tube = [plex_id for plex_id, plex_info in all_plex_info.items() if plex_info['tube'] == i_tube]

        # load existing primers
        #existing_primer = config['existing_primer']
        existing_primer = []

        # randomly select one primer set (one pair each plex)
        curr_index = {}
        curr_fp = {}
        curr_rp = {}

        curr_index_lim = {}     # ada
        curr_badness = 0
        for plex_id in curr_tube:
            # randomly select one primer pair
            i = floor(rand() * len(all_plex_info[plex_id]['optimize']))
            curr_index[plex_id] = i

            curr_index_lim[plex_id] = len(all_plex_info[plex_id]['optimize'])  # ada

            fp = all_plex_info[plex_id]['optimize'][i][0]
            rp = all_plex_info[plex_id]['optimize'][i][1]
            curr_fp[plex_id] = fp['seq']
            curr_rp[plex_id] = rp['seq']
            curr_badness += fp['badness'] + rp['badness']
        inter_badness, comp_badness = PrimerSetBadnessFast(list(curr_fp.values()), list(curr_rp.values()), existing_primer)
        curr_badness += inter_badness
        print('initial loss = %.3f' % curr_badness)

        ## ada

        start_index_seq = index_to_seq(list(curr_index.values()))

        combinations = pd.DataFrame(
            {
                "combination": start_index_seq,
                "badness": curr_badness,
            },
            index=[0],
            #ignore_index=True,
        )
        #measured_index_list.append(start_index_seq)

        ## ada


        # optimization
        learning_curve = []
        SATemp = InitSATemp
        for step in range(NUMSTEPS + ZEROSTEPS):
            cost = 0
            #print('SA temperature = %.3f, SADDLE Loss = %.3f' % (SATemp, curr_badness))
            ## ada
            measured_combination_set = set(combinations["combination"])
            bottom_loss = combinations["badness"].min()

            print('SA temperature = %.3f, SADDLE Loss = %.3f' % (SATemp, bottom_loss))

            bottom_inds = combinations["badness"] <= (1+1/(1+step))* bottom_loss# ;c_threshold * bottom_loss
            #bottom_inds = combinations["badness"] <= c_threshold * bottom_loss
            parents = np.resize(
                combinations["combination"][bottom_inds], #.to_numpy
                batch_size,
            )
            ind_sequences = {}
            #ind_sequences.update(zip(combinations["combination"][bottom_inds], combinations["badness"][bottom_inds]))
            while cost <= cost_budge:
                #parents = _recombine_population_c(parents ,0.1)
                parents = _recombine_population(parents, recomb_rate)
                for i in range(0, len(parents), eval_batch_size):
                    roots_loss = []
                    roots = parents[i : i + eval_batch_size]
                    roots_inds = []
                    for r in roots:
                        if r not in combinations['combination'].tolist():
                            _curr_badness, _, = compute_badness(all_plex_info, curr_tube, r, existing_primer)

                            roots_loss.append(_curr_badness)
                            cost += 1
                        else:
                            roots_loss.append(combinations[combinations['combination']==r]['badness'].values[0])
                            cost += 1
                    nodes = list(enumerate(roots))
                    while (
                        len(nodes) > 0
                        and cost + eval_batch_size < batch_size
                    ):
                        child_idxs = []
                        children = []
                        while len(children) < len(nodes):
                            idx, node_str = nodes[len(children) - 1]
                            node = seq_to_index(node_str)
                            node_amp = dict(zip(curr_tube,node))
                            child = generate_new_combination(node_amp,1/len(node),curr_index_lim)
                            # print(child)
                            child_str = index_to_seq(child)

                            if (child_str not in measured_combination_set
                                    and child_str not in ind_sequences):
                                child_idxs.append(idx)
                                children.append(child_str)
                        children_badness = []
                        for ch in children:
                            if ch not in combinations['combination'].tolist():
                                _curr_badness, _, = compute_badness(all_plex_info, curr_tube, ch, existing_primer)

                                children_badness.append(_curr_badness)
                                cost += 1
                            else:
                                children_badness.append(combinations[combinations['combination']==ch]['badness'].values[0])
                                cost += 1
                            # print("ch", ch)

                        #children_badness =
                        ind_sequences.update(zip(children, children_badness))
                        nodes = []
                        for idx, child, badness in zip(child_idxs, children, children_badness):
                            if badness <= roots_loss[idx]:
                                nodes.append((idx, child))

                # print("get here")
                # print("cost", cost)
                # print("i", i)
                # print("len(nodes):",len(nodes))
                # print("")
            if len(ind_sequences)==0:

                print("debug")
                raise ValueError(
                    "No new combination generated."
                )
            new_combinations = np.array(list(ind_sequences.keys()))
            pred_badness = np.array(list(ind_sequences.values()))
            sorted_order = np.argsort(pred_badness)[:top_num] ############

            # combinations = combinations.append(
            #     pd.DataFrame(
            #         {
            #             "combination": new_combinations[sorted_order],
            #             "badness": pred_badness[sorted_order],
            #         }
            #     )
            # )
            combinations = pd.concat([combinations, pd.DataFrame(
                {
                        "combination": new_combinations[sorted_order],
                        "badness": pred_badness[sorted_order],
                    }
            )], ignore_index=True)
            
        # save optimization result
        for n, (plex_id, i) in enumerate(curr_index.items()):
            all_plex_info[plex_id]['fP'] = all_plex_info[plex_id]['optimize'][i][0]
            all_plex_info[plex_id]['rP'] = all_plex_info[plex_id]['optimize'][i][1]
            all_plex_info[plex_id]['fP_badness'] = comp_badness[0][n]
            all_plex_info[plex_id]['rP_badness'] = comp_badness[1][n]

        # plot learning curve
        fig = plt.figure(figsize=(8, 5))
        plt.plot(learning_curve)
        plt.title('Primer dimer optimization (pool-%d)' % i_tube, fontsize=14)
        plt.xlabel('step', fontsize=14)
        plt.ylabel('SADDLE Loss', fontsize=14)
        save_path = os.path.join(config['out_path'], '%s_pool-%d.png' % (config['title'], i_tube))
        plt.savefig(save_path, bbox_inches='tight')
        plt.close()
        print('Pool %d primer dimer optimization saved as %s' % (i_tube, save_path))

        all_lc.append(learning_curve)
    
    print('primer dimer optimization finished in %.3fs' % (time()-tik))
    # with open('all_plex_info_optimize.pkl', 'wb') as f:
    #     pickle.dump(all_plex_info, f)
    print('Done\n')
    return all_plex_info, all_lc


def save(all_plex_info, config):
    '''
    save to csv
    Input:
        all_plex_info: output of optimize_2(all_plex_info, config), or load from 'all_plex_info_optimize_2.pkl'
        config: design configurations, load from 'config.json'
    Output:
        DataFrame for each pool
    '''
    print('Saving...')
    n_tube = 1
    #n_tube = 2

    l_fP_adp = len(config['fP_prefix'])
    l_rP_adp = len(config['rP_prefix'])

    plex_name = [[] for _ in range(n_tube)] # don't use [[]]*n_tube, this will create linked lists
    id_name = [[] for _ in range(n_tube)]

    # sequence specific part of primers
    fp = [[] for _ in range(n_tube)]
    rp = [[] for _ in range(n_tube)]

    # coordinates
    start = [[] for _ in range(n_tube)]
    insert_start = [[] for _ in range(n_tube)]
    insert_end = [[] for _ in range(n_tube)]
    end = [[] for _ in range(n_tube)]

    # sequences
    amp = [[] for _ in range(n_tube)]
    insert = [[] for _ in range(n_tube)]

    # full primer (with adapters and UMI)
    fp_full = [[] for _ in range(n_tube)]
    rp_full = [[] for _ in range(n_tube)]

    # primer interaction badness
    fp_bad = [[] for _ in range(n_tube)]
    rp_bad = [[] for _ in range(n_tube)]

    # BLAST each primer
    fp_blast = [[] for _ in range(n_tube)]
    rp_blast = [[] for _ in range(n_tube)]

    # amplicon BLAST result
    amplicon_blast = [[] for _ in range(n_tube)]

    # assign plexes to tubes and fill out informations
    #print('amplicons containing SNPs or iSNVs in primers:')
    for plex_id, plex_info in all_plex_info.items():
        i = int(plex_id.split('_')[-1])
        tube = i % n_tube
        if tube == 0:
            tube = n_tube
        plex_name[tube-1].append(plex_id)
        
        curr_fp = plex_info['fP']
        curr_rp = plex_info['rP']

        curr_fp_seq = curr_fp['seq'][l_fP_adp:]
        curr_rp_seq = curr_rp['seq'][l_rP_adp:]
        # if curr_fp_seq != curr_fp_seq.lower() or curr_rp_seq != curr_rp_seq.lower():
        #     print(plex_id)

        fp[tube-1].append(curr_fp['seq'][l_fP_adp:])
        rp[tube-1].append(curr_rp['seq'][l_rP_adp:])

        pc = plex_info['primers_coords'] # primer design coords ([fp_start, fp_end, rp_start, rp_end])
        curr_start = pc[1] - curr_fp['dist'] - curr_fp['primer_len'] + 1
        curr_insert_start = curr_start + curr_fp['primer_len']
        curr_insert_end = pc[2] + curr_rp['dist'] - 1
        curr_end = curr_insert_end + curr_rp['primer_len']
        start[tube-1].append(curr_start)
        insert_start[tube-1].append(curr_insert_start)
        insert_end[tube-1].append(curr_insert_end)
        end[tube-1].append(curr_end)

        curr_context_seq = plex_info['context_seq']
        amp[tube-1].append(curr_context_seq[curr_start-pc[0]:curr_end-pc[0]+1])
        insert[tube-1].append(curr_context_seq[curr_insert_start-pc[0]:curr_insert_end-pc[0]+1])

        fp_full[tube-1].append(curr_fp['seq'])
        rp_full[tube-1].append(curr_rp['seq'])

        fp_bad[tube-1].append(plex_info['fP_badness'])
        rp_bad[tube-1].append(plex_info['rP_badness'])


    # save to file
    df = []
    for n in range(n_tube):
        df.append(pd.DataFrame({
            'amp_id': plex_name[n], 
            'fP': fp[n], 
            'rP': rp[n], 
            'start': start[n], 
            'insert_start': insert_start[n], 
            'insert_end': insert_end[n], 
            'end': end[n], 
            'amp': amp[n], 
            'insert': insert[n], 
            'fP_full': fp_full[n], 
            'rP_full': rp_full[n], 

        }))

    for n, d in enumerate(df):
        save_path = os.path.join(config['out_path'], '%s_pool-%d.csv' % (config['title'], n+1))
        d.to_csv(save_path, index=False)
        print('Primer pool %d saved as %s' % (n+1, save_path))
    return df


def tiling(ref_path, out_path, title, max_amp_len, min_amp_len, w_egc, w_lc, w_ns, w_var, temperature, 
    salinity, dG_max, min_GC, max_GC, min_complexity, max_len, fP_prefix, rP_prefix, seed, threads):
    '''
    Design tiled amplicons. 
    Input:
        ref_path: Path to the Olivar reference file (.olvr).
        out_path: Output directory [./].
        title: Name of this design [olivar-design].
        max_amp_len: Maximum amplicon length [420].
        min_amp_len: Minimum amplicon length. 0.9*{max-amp-len} if None. Minimum 120.
        w_egc: Weight for extreme GC content [1.0].
        w_lc: Weight for low sequence complexity [1.0].
        w_ns: Weight for non-specificity [1.0].
        w_var: Weight for variations [1.0].
        temperature: PCR annealing temperature [60.0].
        salinity: Concentration of monovalent ions in units of molar [0.18].
        dG_max: Maximum free energy change of a primer in kcal/mol [-11.8].
        min_GC: Minimum GC content of a primer [0.2].
        max_GC: Maximum GC content of a primer [0.75].
        min_complexity: Minimum sequence complexity of a primer [0.4].
        max_len: Maximum length of a primer [36].
        fP_prefix: Prefix of forward primer. Empty string '' by default.
        rP_prefix: Prefix of reverse primer. Empty string '' by default.
        seed: Random seed for optimizing primer design regions and primer dimer [10].
        threads: Number of threads.
    '''
    config = {
        'ref_path': ref_path, 
        'out_path': out_path, 
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

    if not os.path.exists(config['out_path']):
        os.makedirs(config['out_path'])

    all_plex_info, risk_arr, gc_arr, comp_arr, hits_arr, var_arr, all_loss, seq_raw = design_context_seq(config)
    all_plex_info_primer = get_primer(all_plex_info, config)
    all_plex_info_optimize, learning_curve = optimize(all_plex_info_primer, config)
    df = save(all_plex_info_optimize, config)

    # save design file
    design_out = {
        'config': config, 
        'risk_arr': risk_arr, 
        'gc_arr': gc_arr, 
        'comp_arr': comp_arr, 
        'hits_arr': hits_arr, 
        'var_arr': var_arr, 
        'all_loss': all_loss, 
        'seq': seq_raw, 
        'df': df, 
        'all_plex_info': all_plex_info, 
        'learning_curve': learning_curve
    }
    design_path = os.path.join(config['out_path'], '%s.olvd' % config['title'])
    with open(design_path, 'wb') as f:
        pickle.dump(design_out, f, protocol=4)
        print('Design file saved as %s' % design_path)


def visualize(design_path):
    '''
    Visualize Olivar designed PCR tiling primers.
    Input:
        design_path: Path to the Olivar design file (.olvd)
    '''
    with open(design_path,  'rb') as f:
        design_out = pickle.load(f)
    config = design_out['config']
    risk_arr = design_out['risk_arr']
    gc_arr = design_out['gc_arr']
    comp_arr = design_out['comp_arr']
    hits_arr = design_out['hits_arr']
    var_arr = design_out['var_arr']
    df = design_out['df']
    print('Configurations:')
    print(json.dumps(config, indent=4))

    r = np.arange(len(risk_arr))

    Plot, Axis = plt.subplots()
    Plot.set_figheight(10)
    Plot.set_figwidth(20)
    plt.subplots_adjust(bottom=0.25) # the bottom space

    # plot risk array
    layer0 = gc_arr
    layer1 = layer0 + comp_arr
    layer2 = layer1 + hits_arr
    layer3 = layer2 + var_arr
    l1 = plt.fill_between(r, layer0, alpha=0.5, color='C0', edgecolor=None, label='extreme GC')
    l2 = plt.fill_between(r, layer1, layer0, alpha=0.5, color='C1', edgecolor=None, label='low complexity')
    l3 = plt.fill_between(r, layer2, layer1, alpha=0.5, color='C2', edgecolor=None, label='non-specificity')
    l4 = plt.fill_between(r, layer3, layer2, alpha=0.5, color='C3', edgecolor=None, label='SNPs')
    l5, = plt.plot([0, 1], [-1, -1], 'k', linewidth=2, label='pool 1 (upper)')
    l6, = plt.plot([0, 1], [-1, -1], 'k', linewidth=2, label='pool 2 (lower)')
    plt.legend(handles=[l1, l2, l3, l4, l5, l6], loc='upper right', fontsize=16)

    #----------- plot primers -----------
    base_offset = (0.5*config['w_egc'] + 0.5*config['w_lc'] + config['w_ns'])/6
    fp_rp_diff = 0.1

    # arrow head
    head_dx = 7
    head_dy = base_offset*7/90

    # pool 1
    fp_offset = 2 * base_offset
    rp_offset = (2-fp_rp_diff) * base_offset
    for i, row in df[0].iterrows():
        fp_start = row['start']-1
        fp_stop = row['insert_start']-1
        rp_start = row['insert_end']
        rp_stop = row['end']
        plt.plot([fp_start, fp_stop], [fp_offset, fp_offset], linewidth=2, color='k')
        plt.plot([fp_stop-head_dx, fp_stop], [fp_offset+head_dy, fp_offset], linewidth=2, color='k')
        plt.plot([rp_start, rp_stop], [rp_offset, rp_offset], linewidth=2, color='k')
        plt.plot([rp_start, rp_start+head_dx], [rp_offset, rp_offset-head_dy], linewidth=2, color='k')

    # pool 2
    fp_offset = base_offset
    rp_offset = (1-fp_rp_diff) * base_offset
    for i, row in df[1].iterrows():
        fp_start = row['start']-1
        fp_stop = row['insert_start']-1
        rp_start = row['insert_end']
        rp_stop = row['end']
        plt.plot([fp_start, fp_stop], [fp_offset, fp_offset], linewidth=2, color='k')
        plt.plot([fp_stop-head_dx, fp_stop], [fp_offset+head_dy, fp_offset], linewidth=2, color='k')
        plt.plot([rp_start, rp_stop], [rp_offset, rp_offset], linewidth=2, color='k')
        plt.plot([rp_start, rp_start+head_dx], [rp_offset, rp_offset-head_dy], linewidth=2, color='k')
    fp_offset += 1
    #----------- plot primers -----------

    plt.title('risk components are stacked together', loc='right', size=16)
    plt.ylim(bottom=0, top=base_offset*6)
    plt.xlabel('position', size=20)
    plt.xticks(fontsize=14)
    plt.ylabel('risk', size=20)
    plt.yticks(fontsize=14)

    slider_color = 'White'

    # Set the axis and slider position in the plot
    axis_position = plt.axes([0.2, 0.1, 0.65, 0.03],
                            facecolor = slider_color)
    slider_position = Slider(axis_position,
                            '', 0.1, len(risk_arr))

    # update() function to change the graph when the
    # slider is in use
    def update(val):
        pos = slider_position.val
        Axis.axis([pos, pos+2000, 0, base_offset*6])
        #Axis.axis([pos, pos+2000, 0, 0.77])
        Plot.canvas.draw_idle()

    # update function called using on_changed() function
    slider_position.on_changed(update)

    # Display the plot
    plt.show()


def validate(primer_pool, BLAST_db, out_path, title, max_amp_len, temperature, threads):
    '''
    Run analysis on a primer pool for multiplex PCR
    Input:
        primer_pool: Path to the csv file of a primer pool. 
            Required columns: "amp_id" (amplicon name), "fP" (sequence of forward primer), "rP" (sequence of reverse primer). 
        BLAST_db: Optional, path to the BLAST database. 
            Note that this path should end with the name of the BLAST database (e.g., "example_input/Human/GRCh38_primary").
        out_path: Output directory [./]. 
        title: Name of validation [olivar-val].
        max_amp_len: Maximum length of predicted non-specific amplicon [1500]. 
        temperature: PCR annealing temperature [60.0]. 
        threads: Number of threads. 
    '''
    if not os.path.exists(out_path):
        os.makedirs(out_path)

    seq_list = []
    seq_names = []
    # df = pd.read_excel(design_path, pool, index_col=None)
    df = pd.read_csv(primer_pool, index_col=False)
    for i, row in df.iterrows():
        seq_list.extend([row['fP'], row['rP']])
        seq_names.extend(['%s_fP' % row['amp_id'], '%s_rP' % row['amp_id']])

    # non-specific simulation
    if BLAST_db:
        df_ns, df_count, all_hits = ns_simulation(BLAST_db, seq_list, seq_names, max_amp_len, threads)
        save_path = os.path.join(out_path, '%s_ns-amp.csv' % title)
        df_ns.to_csv(save_path, index=False)
        print('Non-specific amplicons saved as %s' % save_path)
        save_path = os.path.join(out_path, '%s_ns-pair.csv' % title)
        df_count.to_csv(save_path, index=False)
        print('Non-specific primer pairs saved as %s' % save_path)
    else:
        all_hits = np.zeros(len(seq_list)) - 1
        print('No BLAST database provided, skipped.')

    # calculate Badness
    #all_rep = [seq2arr(p) for p in seq_list]
    _, all_bad = PrimerSetBadnessFast(seq_list)
    all_bad = all_bad[0]

    # calculate dG
    gen = design.primer_generator(temperature=temperature, salinity=0.18)
    all_dG = [gen.dG_init+gen.StacksDG(p) for p in seq_list]

    df_val = pd.DataFrame({
        'name': seq_names, 
        'seq': seq_list, 
        'length': [len(p) for p in seq_list], 
        '%GC': [basic.get_GC(p) for p in seq_list], 
        'complexity': [basic.get_complexity(p) for p in seq_list], 
        'dG (%dC)' % temperature: all_dG, 
        'self_align': [design.WAlignScore(p) for p in seq_list], 
        'dimer_score': all_bad, 
        'BLAST_hits': all_hits
    })
    save_path = os.path.join(out_path, '%s.csv' % title)
    df_val.to_csv(save_path, index=False)
    print('Validation file saved as %s' % save_path)
