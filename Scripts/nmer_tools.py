#!/bin/env python

import itertools as it
import numpy as np
import matplotlib as plt
import fasta_subseq_2 as fa
import itertools as it
import copy as cp

code_fn = {'A':0,'T':1,'G':2,'C':3}
inv_code_fn = ['A','T','G','C']

def build_nmer_profile(pks,win,n,fastadb):
    bases = ['A','T','G','C']
    rc = {'A':'T','C':'G','G':'C','T':'A'}
    nmer_decode_pr = it.product(bases,repeat=n)
    stop = False
    nmer_decode = []
    while not stop:
        try:
            nmer = nmer_decode_pr.next()
            nmer_rc = [rc[x] for x in reversed(nmer)]
            if "".join(nmer_rc) not in nmer_decode:
                nmer_decode.append("".join(nmer))
        except StopIteration:
            stop = True
    nmer_code = dict([(x,i) for (i,x) in enumerate(nmer_decode)])
    nmer_code.update([(x,i) for (i,x) in enumerate(["".join([rc[y] for y in reversed(list(z))]) for z in nmer_decode])])
    profiles = []
    locations = []
    totals = np.zeros(len(nmer_decode))
    for p in pks:
        nmer_ar = np.zeros(len(nmer_decode))
        seq = fastadb[p['chr']][p['summit']-int(win/2):p['summit']+int(win/2)]
            
        for i in range(0,len(seq)-(n-1)):
            subseq = seq[i:i+n]
            if re.search('[^ACGTacgt]',subseq):
                continue
            nmer_ar[nmer_code[subseq]] += 1
        profiles.append(nmer_ar)
        totals += nmer_ar
    return (np.array(profiles), totals,nmer_decode)

# relate nmers to matrices by:
#  1.) calculate matrix score for all nmers at a position (1->n) in matrix
#  2.) take inverse score rank for each nmer 
#  3.) sum for all positions and norm by mtx length for each nmer to give final score

def score_to_ranks(array):
    score_sorted = sorted([(i,s) for (i,s) in enumerate(array)],cmp=(lambda x,y:cmp(x[1],y[1])))
    sorted_order = sorted([(i[0],r,i[1]) for (r,i) in enumerate(score_sorted)],cmp=(lambda x,y:cmp(x[0],y[0])))
    return np.array([x[1] for x in sorted_order])

def count_mtx_to_p_mtx(cm_file,pcounts=True):
    cm_in = open(cm_file)
    freq_mtx = None
    pcount_mtx = None
    for ln in cm_in:
        ln_mtch = re.match("([ATGC])\s+\|((\s\w+)+)",ln)
        base = ln_mtch.group(1)
        ln_counts = np.array(ln_mtch.group(2).split(),dtype=np.float32)
        if freq_mtx == None:
            freq_mtx = np.zeros((len(ln_counts),len(code_fn.keys())))
            pcount_mtx = np.zeros((len(ln_counts),len(code_fn.keys())))
        freq_mtx[:,code_fn[base]] = ln_counts
    if pcounts:
        pcount_mtx = np.sqrt(freq_mtx + 1) / 4
    freq_mtx += pcount_mtx
    freq_mtx /= (freq_mtx.sum() + pcount_mtx.sum())
    freq_mtx = np.log(freq_mtx)
    return freq_mtx

def mtx_nmer_ranks(nmers,nmer_code,freq_mtx,scale=False):
    nmer_mtxs = {}
    nmer_scores = np.zeros(len(nmers))
    nmer_sub_scores = []
    k_max = []
    for nmer in nmers:
        nmer_rc = fa.revcomp(nmer)
        m = np.zeros((len(inv_code_fn),len(nmer)))
        m_rc = np.zeros((len(inv_code_fn),len(nmer)))
        for (i,(b,b_rc)) in enumerate(zip(nmer,nmer_rc)):
            m[code_fn[b],i] = 1
            m_rc[code_fn[b_rc],i] = 1
        nmer_mtxs[nmer] = (m,m_rc)
    
    min_nmer = 5
    max_nmer = len(nmers[0])
    max_score = sum([max(freq_mtx[x,:]) for x in arange(freq_mtx.shape[0])])
    if len(freq_mtx[:,0]) > max_nmer:
        for k in arange(min_nmer,max_nmer+1):
            sub_max = -10000
            k_scores = []
            for z in arange(0,(max_nmer - k)+1):
                for j in arange(0,(len(freq_mtx[:,0]) - k)+1):
                    submtx = freq_mtx[j:j+k,:].T
                    sub_max = max(sub_max,sum([max(submtx[:,x]) for x in arange(submtx.shape[1])]))
                    sub_scores = [(n,max(sum(submtx*m[0][:,z:z+k]),sum(submtx*m[1][:,z:z+k]))) for (n,m) in nmer_mtxs.items()]
                    sub_arr = np.zeros(len(sub_scores))
                    for (n,m) in sub_scores:
                        sub_arr[nmer_code[n]] = m
                    k_scores.append(sub_arr)
                    sub_scores_srt = sorted(sub_scores, cmp=(lambda x,y:cmp(x[1],y[1])))
                    sub_scores_decoded = sorted([(nmer_code[n[0]],n[1],n[0]) for (i,n) in enumerate(sub_scores)],cmp=(lambda x,y:cmp(x[0],y[0])))
                    nmer_scores += np.array([x[1] for x in sub_scores_decoded])
            nmer_sub_scores.extend(k_scores)
            k_arr = np.array(k_scores)
            k_max_scrs = [max(k_arr[:,x]) for x in arange(len(k_scores[0]))]
            k_max_ranks = score_to_ranks(k_max_scrs)
            if scale:
                k_max_ranks *= (sub_max / max_score)
            k_max.append(k_max_ranks)
    else:
        for j in arange(0,len(nmer) - len(freq_mtx[:,0])):
            l = len(freq_mtx[:,0])
            sub_scores = [(n,sum(freq_mtx.T*m[0][j:j+l,:]) + (freq_mtx.T*m[1][j:j+l,:])) for (n,m) in nmer_mtxs.items()]
            sub_arr = np.zeros(len(sub_scores))
            for (n,m) in sub_scores:
                sub_arr[nmer_code[n]] = m
            nmer_sub_scores.append(sub_arr)
            sub_scores_srt = sorted(sub_scores, cmp=(lambda x,y:cmp(x[1],y[1])))
            sub_scores_decoded = sorted([(nmer_code[n[0]],n[1],n[0]) for (i,n) in enumerate(sub_scores)],cmp=(lambda x,y:cmp(x[0],y[0])))
            nmer_scores += np.array([x[1] for x in sub_scores_decoded])

    nmer_sub_scores = np.vstack(nmer_sub_scores)
    nmer_max_subscore = np.array([max(nmer_sub_scores[:,i]) for i in np.arange(len(nmer_sub_scores[0]))])
    nmer_mxss_rank = score_to_ranks(nmer_max_subscore)
    k_max = np.vstack(k_max)
    ss_rank_maxs = [max(k_max[:,x]) for x in arange(k_max.shape[1])]
                    
    nmer_score_srt = sorted([(i,s) for (i,s) in enumerate(nmer_scores)],cmp=(lambda x,y:cmp(x[1],y[1])))
    nmer_ranks = [x[1] for x in sorted([(n[0],i) for (i,n) in enumerate(nmer_score_srt)],cmp=(lambda x,y:cmp(x[0],y[0])))]
    nmer_scores /= freq_mtx.shape[0]
    nmer_scores /= max(abs(nmer_scores))
    return (nmer_ranks,nmer_scores,nmer_mxss_rank,ss_rank_maxs)
    
def nmer_PWM_rank_map(n,matrix_names,*matrix_files):
    pwm_rank_map = {}
    pwm_score_map = {}
    pwm_maxscore_map = {}
    pwm_maxsubrank_map = {}
    bases = ['A','T','G','C']
    rc = {'A':'T','C':'G','G':'C','T':'A'}
    nmer_decode_pr = it.product(bases,repeat=n)
    stop = False
    nmer_decode = []
    while not stop:
        try:
            nmer = nmer_decode_pr.next()
            nmer_rc = [rc[x] for x in reversed(nmer)]
            if "".join(nmer_rc) not in nmer_decode:
                nmer_decode.append("".join(nmer))
        except StopIteration:
            stop = True
    nmer_code = dict([(x,i) for (i,x) in enumerate(nmer_decode)])
    nmer_code.update([(x,i) for (i,x) in enumerate(["".join([rc[y] for y in reversed(list(z))]) for z in nmer_decode])])
    profiles = []
    for (j,m) in enumerate(matrix_files):
        #print (i,m)
        print matrix_names[j]
        (ranks,scores,maxscores,maxsubrank) = mtx_nmer_ranks(nmer_decode,nmer_code,count_mtx_to_p_mtx(m))
        pwm_rank_map[matrix_names[j]] = dict([(nmer_decode[i],r) for (i,r) in enumerate(ranks)])
        pwm_score_map[matrix_names[j]] = dict([(nmer_decode[i],r) for (i,r) in enumerate(scores)])
        pwm_maxscore_map[matrix_names[j]] = dict([(nmer_decode[i],r) for (i,r) in enumerate(maxscores)])
        pwm_maxsubrank_map[matrix_names[j]] = dict([(nmer_decode[i],r) for (i,r) in enumerate(maxsubrank)])
        print map(len,(ranks,scores,maxscores,maxsubrank))
    return (pwm_rank_map,pwm_score_map,pwm_maxscore_map,pwm_maxsubrank_map)

def plot_nmer_slopes(nmer_set,n,mtx_names,mtx_files,ticks=True):
    plt.subplot2grid((5,5),(0,0),colspan=4,rowspan=3)
    counts = np.vstack([smooth(x[3],20) for x in nmer_set]).T
    plt.imshow(counts,cmap='PuRd',aspect='auto',interpolation='none',vmax=0.2)
    plt.subplot2grid((5,5),(3,0),colspan=4,rowspan=1)
    (rank_map,score_map,max_map,sub_max_map) = nmer_PWM_rank_map(n,mtx_names,mtx_files)
    plot_nmer_mtx_ranks(nmer_set,sub_max_map)
    
def plot_nmer_mtx_ranks(nmer_set,rank_map):  
    f_names = []
    nmers = [x[1] for x in nmer_set]
    mr_arr = np.zeros((len(rank_map.keys()),len(nmers)))
    for (i,(n,rmap)) in enumerate(rank_map.items()):
        f_names.append(n)
        print n
        print rmap
        mr_arr[i,:] = np.array([rmap[n] for n in nmers])
        #print mr_arr
        vm=len(rank_map.values()[0].values())
        plt.imshow(mr_arr,cmap='GnBu',aspect='auto',interpolation='none', vmin=vm/2, vmax=vm)
    plt.yticks(arange(len(f_names)),f_names)
    plt.xticks(arange(len(nmer_set)),[x[1] for x in nmer_set],rotation='vertical')

def correlate_nmers(peaks,win,fdr,rand_reps,n,sort_by):
    sorted_rows = sp.plot_by_diff_rank(False,sort_by,False,peaks,ant_w1_norm,pst_w1_norm,wh1_w1_norm)
    sorted_peaks = [x[0] for x in sorted_rows]
    diff_scrs = np.array([x[1] for x in sorted_rows])
    (nmers,profile,decode) = build_nmer_profile(sorted_peaks,win,n,dm_fa)
    significant = []
    all_nmer_corr = []
    fdr_cutoffs = []
    for (i,nmer) in enumerate(decode):
        nmer_counts = nmers[:,i]
        nmer_mtx = np.vstack([nmer_counts,np.ones(len(nmer_counts))]).T
        (sample_m,sample_c) = np.linalg.lstsq(nmer_mtx,diff_scrs)[0]
        linreg_pval = stats.linregress(nmer_counts,diff_scrs)[3]
        if i % 200 == 0:
            print "%d %s sample_m: %f linreg_pval: %f" % (i,nmer,sample_m,linreg_pval)
        all_nmer_corr.append((linreg_pval,nmer,sample_m,nmer_counts))
    all_corr_sorted_linreg = sorted(all_nmer_corr,cmp=(lambda x,y:cmp(y[2],x[2])))
    m = len(all_corr_sorted_linreg)
    significant_linreg = [x for (i,x) in enumerate(all_corr_sorted_linreg) if x[0] < ((float(i)/m)*fdr)]
    return (significant_linreg,all_corr_sorted_linreg,diff_scrs)

def build_nmer_position_profile(pks,win,n,fastadb,loc_win=10):
    bases = ['A','T','G','C']
    rc = {'A':'T','C':'G','G':'C','T':'A'}
    nmer_decode_pr = it.product(bases,repeat=n)
    stop = False
    nmer_decode = []
    while not stop:
        try:
            nmer = nmer_decode_pr.next()
            nmer_rc = [rc[x] for x in reversed(nmer)]
            if "".join(nmer_rc) not in nmer_decode:
                nmer_decode.append("".join(nmer))
        except StopIteration:
            stop = True
    nmer_code = dict([(x,i) for (i,x) in enumerate(nmer_decode)])
    nmer_code.update([(x,i) for (i,x) in enumerate(["".join([rc[y] for y in reversed(list(z))]) for z in nmer_decode])])
    profiles = []
    locations = []
    totals = np.zeros(len(nmer_decode))
    positions = np.zeros((len(nmer_decode),len(pks),win/loc_win),np.int8)
    for (k,p) in enumerate(pks):
        nmer_ar = np.zeros(len(nmer_decode))
        seq = fastadb[p['chr']][p['summit']-int(win/2):p['summit']+int(win/2)]
        posns = []
        if k % 50 == 0:
            print (k,p)
        for i in range(0,len(seq)-(n-1)):
            subseq = seq[i:i+n]
            if re.search('[^ACGTacgt]',subseq):
                continue
            nmer_ar[nmer_code[subseq]] += 1
            positions[nmer_code[subseq],k,int(i/loc_win)] = 1
        profiles.append(nmer_ar)
        totals += nmer_ar
    return (np.array(profiles),totals,nmer_decode,positions)

