import fasta_subseq_2 as fa
import cPickle as pck
import seq_plotmethods as sp
import itertools as it
import tables as tb
import numpy as np
from scipy import stats
import patser_tools as pat
import annotation_track
import matplotlib.pyplot as plt
import random as rnd
import os
import re
from __future__ import division

# tag-normalized tag counts
ant_dpse_h5 = tb.openFile("/Users/cwbrown/Data/Illumina/12_7_2012/h5_ext300/Ant_Gt_ChIP_e300_dpse.h5", mode='r+')
pst_dpse_h5 = tb.openFile("/Users/cwbrown/Data/Illumina/12_7_2012/h5_ext300/Post_Gt_ChIP_e300_dpse.h5", mode='r+')
cmb_dpse_h5 = tb.openFile("/Users/cwbrown/Data/Illumina/12_7_2012/h5_ext300/Combo_Gt_ChIP_e300_dpse.h5", mode='r+')
wh1_dpse_h5 = tb.openFile("/Users/cwbrown/Data/Illumina/12_7_2012/h5_ext300/Whole1_Gt_ChIP_e300_dpse.h5", mode='r+')
wh2_dpse_h5 = tb.openFile("/Users/cwbrown/Data/Illumina/12_7_2012/h5_ext300/Whole2_Gt_ChIP_e300_dpse.h5", mode='r+')
wh1_h5 = tb.openFile("/Users/cwbrown/Data/Illumina/12_7_2012/h5_ext300/Whole1_Gt_ChIP_dmel_e300.h5", mode='r+')

# make chipExpt objects for all above (tag-normed and raw counts)
ant_w1_norm = sp.chipExpt(ant_dpse_h5.root.normalized_tags, name="ant_cn",chr_ext="norm")
pst_w1_norm = sp.chipExpt(pst_dpse_h5.root.normalized_tags, name="pst_cn",chr_ext="norm")
cmb_w1_norm = sp.chipExpt(cmb_dpse_h5.root.normalized_tags, name="cmb_cn",chr_ext="norm")
wh2_w1_norm = sp.chipExpt(wh2_dpse_h5.root.normalized_tags, name="wh2_cn",chr_ext="norm")
wh1_w1_norm = sp.chipExpt(wh1_h5.root.tag_counts.exttag_normed_tags,"Whole Rep1 Gt ChIP, tag-normalized")

pks_in = open("/Users/cwbrown/Data/slicing_paper_analysis/Fig5_sequence_correlates/peak_sets.pck")
(ant_fdr1_pks,pst_fdr1_pks,fdr40_ap_nodiff_pks,byexp_top1000_un50_pks) = pck.load(pks_in)
pks_in.close()

rsq_out = open("/Users/cwbrown/Data/slicing_paper_analysis/Fig5_sequence_correlates/nmer_rsq.pck","w")
abs_AP_diff_nmers = correlate_nmers(byexp_top1000_un50_pks,500,.1,30000,6,"abs")
abs_v_AP_diff_nmers = correlate_nmers(byexp_top1000_un50_pks,500,.1,30000,6,"abs_v")
#r_AP_diff_nmers = correlate_nmers(byexp_top1000_un50_pks,500,.1,50000,6,"rel")
nmer_rsq = (abs_AP_diff_nmers,abs_v_AP_diff_nmers)#,r_AP_diff_nmers)
pck.dump(nmer_rsq,rsq_out)
rsq_out.close()

def correlate_nmers(peaks,win,fdr,rand_reps,n,sort_by):
    sorted_rows = sp.plot_by_diff_rank(False,sort_by,False,peaks,ant_w1_norm,pst_w1_norm,wh1_w1_norm)
    sorted_peaks = [x[0] for x in sorted_rows]
    diff_scrs = np.array([x[1] for x in sorted_rows])
    (nmers,profile,decode) = build_nmer_profile(sorted_peaks,win,6,dm_fa)
    print len(decode)
    significant = []
    all_nmer_corr = []
    fdr_cutoffs = []
    rand_sets_mtx = np.zeros((len(diff_scrs),rand_reps))
    rand_set_linfit = np.zeros((len(decode),rand_reps))
    for x in range(0,rand_reps):
        if (x % 1000) == 0:
            print "%d random sets" % x
        rand_sets_mtx[:,x] = np.random.permutation(diff_scrs)
        
    for (i,nmer) in enumerate(decode):
        nmer_counts = nmers[:,i]
        #corr = np.corrcoef(np.array([nmer_counts,rset]))[0,1]
        #rand_set_corr[i].append(corr)
        nmer_mtx = np.vstack([nmer_counts,np.ones(len(nmer_counts))]).T
        linfit = np.linalg.lstsq(nmer_mtx,rand_sets_mtx)
        coeffs = linfit[0][0]
        rand_set_linfit[i,:] = coeffs
        #fdr_ucut = np.percentile(coeffs,100-(fdr/2.0))
        #fdr_lcut = np.percentile(coeffs,fdr/2.0)
        q1 = np.percentile(coeffs,25)
        q3 = np.percentile(coeffs,75)
        #plt.subplot(5,1,i)
        #plt.hist(coeffs,bins=100)
        #plt.axvline(fdr_cut)
        #plt.axvline(-fdr_cut)
        (sample_m,sample_c) = np.linalg.lstsq(nmer_mtx,diff_scrs)[0]
        #l_pval = (len(np.where(coeffs <= sample_m)[0])/len(coeffs))*2 # pval = density of emp dist below sample m, x2 since two-tailed test
        #h_pval = (len(np.where(coeffs >= sample_m)[0])/len(coeffs))*2 # pval = density of emp dist below sample m, x2 since two-tailed test
        #empirical_pval = min((l_pval,h_pval))
        (ntest_stat,ntest_pval) = stats.normaltest(coeffs)
        gkde=stats.gaussian_kde(coeffs)
        l_cont_pval = gkde.integrate_box_1d(-inf,sample_m) * 2
        h_cont_pval = gkde.integrate_box_1d(sample_m,inf) * 2
        kde_pval = min([l_cont_pval,h_cont_pval])
        #if ((i % 10 == 0) in range(5)) and (i <=50):
        #    gkde=stats.gaussian_kde(coeffs)
        #    plt.subplot(5,1,(i/10))
        #    plt.hist(coeffs,bins=50,normed=True)
        #    s_pts = np.linspace(min(coeffs)-1,max(coeffs)+1)
        #    s_pts_pdf = gkde.evaluate(s_pts)
        #    plt.plot(s_pts,s_pts_pdf,"r-",label='Gaussian KDE')
        #    plt.plot(s_pts,stats.norm.pdf(s_pts),'g-',label="Normal")
        #    plt.axvline(sample_m)
        #    plt.legend()
        #if (sample_m > fdr_ucut) or (sample_m < fdr_lcut):
        #    print "SIGNIFICANT NMER: id:%d seq:%s Rsq:%f ucut:%f lcut:%f" % (i,nmer,sample_m,fdr_ucut,fdr_lcut)
        #    significant.append((i,nmer,sample_m,sample_c,nmer_counts,diff_scrs))
        #if i % 50 == 0:
        print "%d %s sample_m: %f kde_pval: %f normaltest_stat: %f normaltest pvalue: %f" % (i,nmer,sample_m,kde_pval,ntest_stat,ntest_pval)
        all_nmer_corr.append((kde_pval,ntest_pval,nmer,sample_m,nmer_counts))
    #all_corr_sorted_emp = sorted(all_nmer_corr,cmp=(lambda x,y:cmp(y[0],x[0])))
    all_corr_sorted_kde = sorted(all_nmer_corr,cmp=(lambda x,y:cmp(y[1],x[1])))
    m = len(all_corr_sorted_kde)
    #significant_emp = [x for (i,x) in enumerate(all_corr_sorted_emp) if x[0] < ((float(i)/m)*fdr)]
    significant_kde = [x for (i,x) in enumerate(all_corr_sorted_kde) if x[0] < ((float(i)/m)*fdr)]
    #all_nmer_cnt_sorted = np.array([x[1] for x in all_corr_sorted]).T
    #plt.imshow(all_nmer_cnt_sorted,cmap='Greys',aspect='auto',interpolation='none') 
    #return ((significant_emp,significant_kde),(all_corr_sorted_emp,all_corr_sorted_kde),diff_scrs)
    return (significant_kde,all_corr_sorted_kde,diff_scrs)

def build_nmer_profile(pks,win,n,fastadb):
    bases = ['A','T','G','C']
    rc = {'A':'T','C':'G','G':'C','T':'A'}
    nmer_decode_pr = it.product(bases,repeat=6)
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
    #print len(nmer_decode)
    nmer_code = dict([(x,i) for (i,x) in enumerate(nmer_decode)])
    print len(nmer_code)
    nmer_code.update([(x,i) for (i,x) in enumerate(["".join([rc[y] for y in reversed(list(z))]) for z in nmer_decode])])
    print len(nmer_code)
    #print nmer_code
    #print nmer_code['CCCCCA']
    #print nmer_decode[4095]
    print nmer_code
    profiles = []
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
