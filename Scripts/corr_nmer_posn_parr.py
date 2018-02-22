#!/usr/bin/env python
#
#
# parallelization of correlate nmer positions (for set of nmers + set of peaks, find correlations between positions of nmers within window win)
#
#

import sys
import os
import re
import time
import numpy as np
import scipy.stats as stats
import multiprocessing as mp
import pickle as pck
import nmer_tools as nm

fasta_file_loc = "/Users/cwbrown/Data/genomes/dmel.fa"

pos_prf = None
p_mtx = None
pval_mtx = None

def _main(args):
    if len(args) != 4:
        print "usage: corr_nmer_posn_parr.py <peaks_pickle> <nmer_length> <pk_window> <corr_window>"
        sys.exit(0)

    pks_in = open(args[0])
    pks = pck.load(pks_in)
    win = args[2]
    n = args[1]
    loc_win = args[3]

    fastadb = fa.FastaDB(fasta_file_loc)

    print "building nmer location matrix...",
    (profiles,totals,nmer_decode,position_profile) = build_nmer_position_profile(pks,win,n,fastadb,loc_win)
    pos_prf = position_profile
    print "done!"
    pval_mtx = np.zeros((position_profile.shape[0],position_profile.shape[0]))
    pos_mtx_tot = position_profile.shape[1] * position_profile.shape[2]
    print "building p_mtx...",
    p_mtx = (np.sum(np.sum(poscr,axis=1),axis=1)/pos_mtx_tot)
    print "done!"
    
    print "building correlation matrix...",
    run_single_nmer = (lambda x: sing_nmer_corr(x,position_profile,p_mtx,overlap_binpvals))
    p = mp.Pool()
    p.map(run_single_nmer,arange(0,position_profile.shape[0]))
    corr_mtx_out = open("nmer_corr_binpval_mtx.pck","w")
    pck.dump(overlap_binpvals,corr_mtx_out)
    corr_mtx_out.close()
    print "All Done, Buddy!"

def sing_nmer_corr(i,pos_prf,p_mtx,pvals_mtx):
    t_st = time.clock()
    m_i = poscr[i,:,:]
    p_i = p_mtx[i]
    m_j = poscr[0:i,:,:]
    m_j_pmtx = p_mtx[0:i] * p_i        
    sum_mtx = m_j + m_i
    ol_mtx = np.array(np.where(sum_mtx > 1))
    overlap_binpvals[i,:i] = [stats.binom_test(len(np.where(ol_mtx[0,:] == k)),pos_mtx_tot,p=m_j_pmtx[k]) for k in np.arange(m_j.shape[0])]
    if i % 20 == 0:
        print "nmer #%d done!  Time elapsed:%f" % (i,(time.clock()-t_st))

if __name__ == "__main__":
    _main(sys.argv[1:]
