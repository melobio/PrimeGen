#!/usr/bin/env python3
# -*- coding: utf-8 -*-


'''
Some basic tools for processing DNA sequences and coordinates
'''


__author__ = 'Michael X. Wang'


import os
CURRENT_DIR = os.path.dirname(__file__)


from collections import Counter

import numpy as np


def revcomp(seq):
    '''
    give the reverse complement of input sequence
    base & number conversion:
        {'A':0, 'T':1, 'C':2, 'G':3}
    input:
        string or array sequence
    output:
        reverse complement of input
    '''
    if isinstance(seq, str):
        complement = {'A':'T', 'T':'A', 'C':'G', 'G':'C', 'N':'N', 
        'a':'t', 't':'a', 'c':'g', 'g':'c', 'n':'n', 
        '0':'1', '1':'0', '2':'3', '3':'2', '4':'4'}
        bases = [complement[base] for base in seq]
        bases.reverse()
        return ''.join(bases)
    elif isinstance(seq, list):
        complement = {'A':'T', 'T':'A', 'C':'G', 'G':'C', 'N':'N', 
        'a':'t', 't':'a', 'c':'g', 'g':'c', 'n':'n', 
        '0':'1', '1':'0', '2':'3', '3':'2', '4':'4', 
        0:1, 1:0, 2:3, 3:2, 4:4}
        bases = [complement[base] for base in seq]
        bases.reverse()
        return bases
    else:
        return -1


def get_GC(seq):
    '''
    get the GC ratio of input sequence
    base & number conversion:
        {'A':0, 'T':1, 'C':2, 'G':3}
    input:
        string or array sequence
    output:
        GC ratio of input sequence
    '''
    return len([base for base in list(seq) if base in ['G','g',3,'C','c',2]])/len(seq)

def get_GC_loss(lower, upper, seq):
    loss = 1
    '''
    get the GC ratio of input sequence
    base & number conversion:
        {'A':0, 'T':1, 'C':2, 'G':3}
    input:
        string or array sequence
    output:
        GC ratio of input sequence
    '''
    ratio = len([base for base in list(seq) if base in ['G','g',3,'C','c',2]])/len(seq)
    if ratio > upper or ratio < lower:
        loss = 1
    else:
        loss = 0
    return loss

def get_GC_loss_2(lower, upper, seq):
    loss = 1
    '''
    get the GC ratio of input sequence
    base & number conversion:
        {'A':0, 'T':1, 'C':2, 'G':3}
    input:
        string or array sequence
    output:
        GC ratio of input sequence
    '''
    ratio = len([base for base in list(seq) if base in ['G','g',3,'C','c',2]])/len(seq)
    if ratio <= upper and ratio >=lower:
        loss = 0
    elif ratio <lower:
        loss = lower - ratio
    elif ratio >upper:
        loss = ratio - upper
    return loss

def get_ATGC_content_loss(lower,upper, seq):
    loss = 1
    '''
    get the GC ratio of input sequence
    base & number conversion:
        {'A':0, 'T':1, 'C':2, 'G':3}
    input:
        string or array sequence
    output:
        GC ratio of input sequence
    '''
    a_content = len([base for base in list(seq) if base in ['A','a',0]])/len(seq)
    t_context = len([base for base in list(seq) if base in ['T','t',2]])/len(seq)
    c_context = len([base for base in list(seq) if base in ['C','c',2]])/len(seq)
    g_context = len([base for base in list(seq) if base in ['G','g',3]])/len(seq)
    if a_content>upper or t_context>upper or c_context>upper or g_context>upper:
        loss = 1
    elif a_content<lower or t_context<lower or c_context<lower or g_context<lower:
        loss = 1
    else:
        loss = 0
    return loss

def get_ATGC_content_loss_2(lower,upper, seq):
    #loss = 0
    '''
    get the GC ratio of input sequence
    base & number conversion:
        {'A':0, 'T':1, 'C':2, 'G':3}
    input:
        string or array sequence
    output:
        GC ratio of input sequence
    '''
    a_content = len([base for base in list(seq) if base in ['A','a',0]])/len(seq)
    t_context = len([base for base in list(seq) if base in ['T','t',2]])/len(seq)
    c_context = len([base for base in list(seq) if base in ['C','c',2]])/len(seq)
    g_context = len([base for base in list(seq) if base in ['G','g',3]])/len(seq)
    #loss a
    if a_content<=upper and a_content >=lower:
        loss_a = 0
    elif a_content < lower:
        loss_a = lower-a_content
    elif a_content > upper:
        loss_a = a_content - upper
    #loss t
    if t_context<=upper and t_context >=lower:
        loss_t = 0
    elif t_context < lower:
        loss_t = lower-t_context
    elif t_context > upper:
        loss_t = t_context - upper
    #loss c
    if c_context<=upper and c_context >=lower:
        loss_c = 0
    elif c_context < lower:
        loss_c = lower-c_context
    elif c_context > upper:
        loss_c = c_context - upper
    #loss g
    if g_context<=upper and g_context >=lower:
        loss_g = 0
    elif g_context < lower:
        loss_g = lower-g_context
    elif g_context > upper:
        loss_g = g_context - upper

    loss = loss_a + loss_t + loss_c + loss_g
    return loss

def get_ATGC_content_count(lower,upper, seq):
    count = 0
    '''
    get the GC ratio of input sequence
    base & number conversion:
        {'A':0, 'T':1, 'C':2, 'G':3}
    input:
        string or array sequence
    output:
        GC ratio of input sequence
    '''
    a_content = len([base for base in list(seq) if base in ['A','a',0]])/len(seq)
    t_context = len([base for base in list(seq) if base in ['T','t',2]])/len(seq)
    c_context = len([base for base in list(seq) if base in ['C','c',2]])/len(seq)
    g_context = len([base for base in list(seq) if base in ['G','g',3]])/len(seq)
    if a_content>upper or a_content<lower:
        count+=1
    
    if t_context>upper or t_context<lower:
        count +=1
    if c_context>upper or c_context<lower:
        count +=1
    if g_context>upper or g_context<lower:
        count += 1
    
    
    return count


def seq2arr(seq_str):
    '''
    convert string sequence to array (list of numbers)
    base & number conversion:
        {'A':0, 'T':1, 'C':2, 'G':3}
    input: 
        string sequence
    output:
        array sequence
    '''
    dic = {'A':0, 'T':1, 'C':2, 'G':3, 'a':0, 't':1, 'c':2, 'g':3}
    seq_list = [dic[s] for s in list(seq_str)]
    return seq_list


def arr2seq(seq_array):
    '''
    convert array to string sequence
    base & number conversion:
        {'A':0, 'T':1, 'C':2, 'G':3}
    input: 
        array sequence
    output:
        string sequence
    '''
    dic = {0:'A', 1:'T', 2:'C', 3:'G'}
    seq_list = [dic[s] for s in seq_array]
    return ''.join(seq_list)


def seq2num(seq):
    '''
    use quaternary to convert array or string sequence to a number
    CAUTION: diffent sequence may output the same number (e.g. 'AACT' and 'AAACT')
    base & number conversion:
        {'A':0, 'T':1, 'C':2, 'G':3}
    input:
        string or array sequence
    output:
        an integer
    '''
    bases = list(seq)
    dic = {'A':0, 'T':1, 'C':2, 'G':3, 
    'a':0, 't':1, 'c':2, 'g':3, 
    0:0, 1:1, 2:2, 3:3}
    bases = [dic[base] for base in bases]
    len_bases = len(seq)
    out_num = 0
    for i in range(len_bases):
        out_num = out_num + 4**(len_bases-i-1)*bases[i]
    return out_num


def num2seq(num, l=None):
    '''
    The reverse function of seq2num
    CAUTION: a single number represents multiple sequences if l=None
    input:
        num: int
        l: length of output sequence; if none, 'A' at sequence beginning is ignored
    output:
        sequence
    '''
    base = ['A', 'T', 'C', 'G']
    seq = []
    if l is None:
        while True:
            new_num = num // 4
            b = num % 4
            seq.append(base[b])
            if new_num == 0:
                break
            num = new_num
        return ''.join(seq[::-1])
    elif num < 4**l:
        for _ in range(l):
            new_num = num // 4
            b = num % 4
            seq.append(base[b])
            num = new_num
        return ''.join(seq[::-1])
    else:
        print('Input number out of range according to input sequence length.')
        return None


def seq2complex(seq_str):
    '''
    Convert string sequence to a list of complex
    '''
    dic = {'A':1, 'T':-1, 'C':1j, 'G':-1j, 'N': 0, 'a':1, 't':-1, 'c':1j, 'g':-1j, 'n': 0}
    seq_list = [dic[s] for s in list(seq_str)]
    return seq_list


def complex2seq(seq_complex):
    '''
    Convert a list of complex back to string
    '''
    seq_str = []
    for c in seq_complex:
        if c.real > c.imag:
            if -c.real < c.imag:
                seq_str.append('A')
            elif -c.real > c.imag:
                seq_str.append('G')
            else:
                seq_str.append('N')
        elif c.real < c.imag:
            if -c.real < c.imag:
                seq_str.append('C')
            elif -c.real > c.imag:
                seq_str.append('T')
            else:
                seq_str.append('N')
        else:
            seq_str.append('N')
    return ''.join(seq_str)


def randseq(l, seed=None):
    '''
    Generate random sequence of length l
    input:
        l: sequence length
        seed: random seed for numpy
    output:
        random sequence of length l
    '''
    if seed is not None:
        np.random.seed(seed)
    return arr2seq(np.random.randint(4, size=l).tolist())


def randarr(l, seed=None):
    '''
    Generate random array of length l
    input:
        l: array length
        seed: random seed for numpy
    output:
        random array of length l
    '''
    if seed is not None:
        np.random.seed(seed)
    return np.random.randint(4, size=l).tolist()


def get_complexity(seq):
    # calculate Shannon entropy for word size of 1, 2 and 3
    hash1 = Counter()
    hash2 = Counter()
    hash3 = Counter()

    # count occurance
    for j in range(len(seq)-2):
        hash1[seq[j]] += 1
        hash2[seq[j:j+2]] += 1
        hash3[seq[j:j+3]] += 1
    hash1[seq[-2]] += 1
    hash1[seq[-1]] += 1
    hash2[seq[-2:]] += 1

    # calculate probability mass function
    p1 = np.array(list(hash1.values()))/(len(seq))
    p2 = np.array(list(hash2.values()))/(len(seq)-1)
    p3 = np.array(list(hash3.values()))/(len(seq)-2)

    SE1 = -np.sum(p1 * np.log2(p1))
    SE2 = -np.sum(p2 * np.log2(p2))
    SE3 = -np.sum(p3 * np.log2(p3))
    return min([SE1/2, SE2/4, SE3/6])
