#!/usr/bin/env python
import sys
import os
import numpy as np
from CODON_TABLE import CODON_TABLE as codons

def translate(seq):
    translated = "".join([codons[seq[3*i:3*i+3]]['1letter'] for i in range(0,len(seq)/3)])
    return translated
   
def get_codons(seq):
    cdns = [seq[i:i+3] for i in range(0,len(seq),3)]
    return cdns

def get_codon_variants(codon,from_nucl=None,to_nucl=None):
    codon_target = codon.upper()
    
    if to_nucl:
        bases = [to_nucl,]
    else:
        bases = ['A','T','G','C']

    p1_block = ["".join([n,codon_target[1],codon_target[2]]) for n in bases if not n == codon_target[0]]
    p2_block = ["".join([codon_target[0],n,codon_target[2]]) for n in bases if not n == codon_target[1]]
    p3_block = ["".join([codon_target[0],codon_target[1],n]) for n in bases if not n == codon_target[2]]

    if from_nucl:
        if from_nucl != codon_target[0]:
            p1_block = []
        if from_nucl != codon_target[1]:
            p2_block = []
        if from_nucl != codon_target[2]:
            p3_block = []

    by_posn = {'p1':[(c,codons[c]) for c in p1_block],\
               'p2':[(c,codons[c]) for c in p2_block],\
               'p3':[(c,codons[c]) for c in p3_block]}
    full_codon_list = p1_block + p2_block + p3_block
    translated_list = [(c,codons[c]) for c in full_codon_list]
    return translated_list,by_posn

def get_synonymous_variants(codon,from_nucl=None,to_nucl=None):
    aa_1letter = codons[codon]['1letter']
    all_variants = get_codon_variants(codon,from_nucl,to_nucl)[0]
    syn_variants = [c for c in all_variants if c[1]['1letter'] == aa_1letter]
    return syn_variants

def get_nonsynonymous_variants(codon,from_nucl=None,to_nucl=None):
    aa_1letter = codons[codon]['1letter']
    all_variants = get_codon_variants(codon,from_nucl,to_nucl)[0]
    nonsyn_variants = [c for c in all_variants if (c[1]['1letter'] != aa_1letter) and (c[1]['1letter'] != "*")]
    return nonsyn_variants

def get_nonsense_variants(codon,from_nucl=None,to_nucl=None):
    aa_1letter = codons[codon]['1letter']
    all_variants = get_codon_variants(codon,from_nucl,to_nucl)[0]
    nonsense_variants = [c for c in all_variants if c[1]['1letter'] == "*"]
    return nonsense_variants

def get_single_base_variants(sequence,codon_target):
    codon_target = codon_target.upper()
    bases = ['A','T','G','C']
    p1_block = ["".join([n,codon_target[1],codon_target[2]]) for n in bases if not n == codon_target[0]]
    p2_block = ["".join([codon_target[0],n,codon_target[2]]) for n in bases if not n == codon_target[1]]
    p3_block = ["".join([codon_target[0],codon_target[1],n]) for n in bases if not n == codon_target[2]]

    blocks = {'p1':p1_block,'p2':p2_block,'p3':p3_block}
    codons_by_var_pos = {'p1':[],'p2':[],'p3':[]}
    seq_codons = get_codons(sequence)
    for i,c in enumerate(seq_codons):
        if c in p1_block:
            codons_by_var_pos['p1'].append((i,c,codons[c]))
        if c in p2_block:
            codons_by_var_pos['p2'].append((i,c,codons[c]))
        if c in p3_block:
            codons_by_var_pos['p3'].append((i,c,codons[c]))

    sorted_by_res = sorted(codons_by_var_pos['p1'] + codons_by_var_pos['p2'] + codons_by_var_pos['p3'],cmp=(lambda x,y:cmp(x[0],y[0])))
    return (blocks,codons_by_var_pos,sorted_by_res)

def print_single_base_variants(sequence,codon_target,line_len=30):
    cdns = get_codons(sequence)
    codon_target = codon_target.upper()
    (blocks,cdns_var_pos,cdns_sorted) = get_single_base_variants(sequence,codon_target)

    top_lines = [[],[],[]]
    bottom_lines = [[],[],[]]
    top_n = 0
    bottom_n = 0

    for (i,mcdn) in enumerate(cdns_sorted):
        if int(i) / 2 == float(i) / 2:
            if mcdn[1] in blocks['p1']:
                p1_not_cdn = [c for c in blocks['p1'] if c != mcdn[1]]
                for l in bottom_lines:
                    l.extend([' ',' ',' '] * (mcdn[0]-bottom_n))
                bottom_lines[2].extend(['*',' ',' ',' ',' ',' '])
                bottom_lines[1].extend(list(p1_not_cdn[0]) + [':',codons[p1_not_cdn[0]]['1letter']," "])
                bottom_lines[0].extend(list(p1_not_cdn[1]) + [':',codons[p1_not_cdn[1]]['1letter']," "])
            elif mcdn[1] in blocks['p2']:
                p2_not_cdn = [c for c in blocks['p2'] if c != mcdn[1]]
                for l in bottom_lines:
                    l.extend([' ',' ',' '] * (mcdn[0]-bottom_n))
                bottom_lines[2].extend([' ','*',' ',' ',' ',' '])
                bottom_lines[1].extend(list(p2_not_cdn[0]) + [':',codons[p2_not_cdn[0]]['1letter']," "])
                bottom_lines[0].extend(list(p2_not_cdn[1]) + [':',codons[p2_not_cdn[0]]['1letter']," "])
            elif mcdn[1] in blocks['p3']:
                p3_not_cdn = [c for c in blocks['p3'] if c != mcdn[1]]
                for l in bottom_lines:
                    l.extend([' ',' ',' '] * (mcdn[0]-bottom_n))
                bottom_lines[2].extend([' ',' ','*',' ',' ',' '])
                bottom_lines[1].extend(list(p3_not_cdn[0]) + [':',codons[p3_not_cdn[0]]['1letter']," "])
                bottom_lines[0].extend(list(p3_not_cdn[1]) + [':',codons[p3_not_cdn[1]]['1letter']," "])
            bottom_n = mcdn[0] + 2
        #print ["".join(p) for p in bottom_lines]

        else:
            if mcdn[1] in blocks['p1']:
                p1_not_cdn = [c for c in blocks['p1'] if c != mcdn[1]]
                for l in top_lines:
                    l.extend([' ',' ',' '] * (mcdn[0]-top_n))
                top_lines[0].extend(['*',' ',' ',' ',' ',' '])
                top_lines[1].extend(list(p1_not_cdn[0]) + [':',codons[p1_not_cdn[0]]['1letter']," "])
                top_lines[2].extend(list(p1_not_cdn[1]) + [':',codons[p1_not_cdn[1]]['1letter']," "])
            elif mcdn[1] in blocks['p2']:
                p2_not_cdn = [c for c in blocks['p2'] if c != mcdn[1]]
                for l in top_lines:
                    l.extend([' ',' ',' '] * (mcdn[0]-top_n))
                top_lines[0].extend([' ','*',' ',' ',' ',' '])
                top_lines[1].extend(list(p2_not_cdn[0]) + [':',codons[p2_not_cdn[0]]['1letter']," "])
                top_lines[2].extend(list(p2_not_cdn[1]) + [':',codons[p2_not_cdn[1]]['1letter']," "])
            elif mcdn[1] in blocks['p3']:
                p3_not_cdn = [c for c in blocks['p3'] if c != mcdn[1]]
                for l in top_lines:
                    l.extend([' ',' ',' '] * (mcdn[0]-top_n))
                top_lines[0].extend([' ',' ','*',' ',' ',' '])
                top_lines[1].extend(list(p3_not_cdn[0]) + [':',codons[p3_not_cdn[0]]['1letter']," "])
                top_lines[2].extend(list(p3_not_cdn[1]) + [':',codons[p3_not_cdn[1]]['1letter']," "])
            top_n = mcdn[0] + 2
            
    aa_line = "".join(np.concatenate([[codons[c]['1letter'],' ',' '] for c in cdns]))
    pos_line = "".join(np.concatenate([list(str(i)) + [' ']*(30-len(list(str(i)))) for (i,c) in enumerate(cdns) if i % 10 == 0]))

    print_strs = []
 
    for i in range(0,((len(top_lines[0])/3)/line_len + 1)):
        st = 3*(i*line_len)
        end = 3*((i+1)*line_len)
        if end > len(top_lines[0]):
            end = len(top_lines[0])
            
        for l in top_lines:
            lst = "".join(l)
            print_strs.append(lst[st:end])
        print_strs.append(aa_line[st:end])
        print_strs.append(sequence[st:end])
        print_strs.append(pos_line[st:end])
        for l in bottom_lines:
            lst = "".join(l)
            print_strs.append(lst[st:end])
        print_strs.append("")
    print_ln = "\n".join(print_strs)
    return print_ln
