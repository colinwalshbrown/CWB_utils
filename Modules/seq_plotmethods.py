#!/usr/bin/env python
#
# seq_plotmethods.py - module of graphical and utility methods for dealing with NGS (ChIP-seq and RNA-seq) data
#

from __future__ import division
import alignment as al
import tables as tb
import matplotlib.pyplot as plt
import itertools as it
import tempfile as tmp
import subprocess as sub
import numpy as np
import re, os
import copy as cp

DMEL_CHRS = ['2L','2R','3L','3R','X','YHet','4','2LHet','2RHet','3LHet','3RHet','XHet','U','Uextra']
PLOT_COLORS = ['r','g','b','c','m','y','k']
PLOT_SHAPES = ['.','o','v','^','<','>','s','+','x','D','d']
LINE_STYLES = ['-','--','-.',':']
BLAST_EXE = "/Users/cwbrown/src/blast3/blastn"
DMEL_GENOME = "/Users/cwbrown/Data/genomes/dmel.fa"

def get_plot_marker():
    for s in PLOT_SHAPES:
        for c in PLOT_COLORS:
            #print s + c
            yield s + c

def get_line_style():
    for s in LINE_STYLES:
        for c in PLOT_COLORS:
            yield c + s

def get_plot_color():
    for p in PLOT_COLORS:
        yield p

def read_macs(chip_hf5,macs_file,expt_name="NoName"):

    try:
        macs_group = chip_hf5.createGroup("/",'macs_peaks','Peaks Called using MACS')
    except tb.exceptions.NodeError:
        print "WARNING: group \'macs_peaks\' exists!"
        macs_group = chip_hf5.root.macs_peaks
    try:
        macs_table = chip_hf5.createTable(macs_group,expt_name + '_peaks',MacsPeak,expt_name + ' Peaks Called using MACS')
    except tb.exceptions.NodeError:
        chip_hf5.removeNode(macs_group,expt_name + '_peaks')
        print "WARNING: table \'"+expt_name+"_peaks\' exists!"
        macs_table = chip_hf5.createTable(macs_group,expt_name + '_peaks',MacsPeak,expt_name + ' Peaks Called using MACS')
        
    for raw_ln in open(macs_file):
        ln = raw_ln.strip().split()
        if len(ln) <= 1 or ln[0][0] == "#":
            #print raw_ln
            continue
        elif ln[0] == 'chr':
            continue
        r = macs_table.row
        r['chr'] = ln[0]
        r['start'] = int(ln[1])
        r['end'] = int(ln[2])
        r['length'] = int(ln[3])
        r['summit'] = int(ln[4])
        r['tags'] = float(ln[5])
        r['log_pval'] = float(ln[6])
        r['fold_enrichment'] = float(ln[7])
        r['log_qval'] = float(ln[8])
        r.append()

    macs_table.flush()
    return macs_table
    
def macs_overlap(print_files,*macs_tables):
    
    """
    DEPRECATED - use overlap() instead
    recursively check for overlapping regions among sets of macs data (read into hf5 table using read_macs method above) 
    """

#    print macs_tables

    if len(macs_tables) < 2:
        return

    overlap_peaks = []

    for r in macs_tables[0].iterrows():
        overlap = []
        ov_peak = None
        a_ov = True
        for t in macs_tables[1:]:
            t_ov = False
            for r2 in t.iterrows():
                if((r['chr'] == r2['chr']) and
                   (((r['start'] <= r2['end']) and (r['start'] >= r2['start'])) or
                    ((r['end'] <= r2['end']) and (r['end'] >= r2['start'])))):
                    #print "r1 => %s[%d:%d]" % (r['chr'],r['start'],r['end'])
                    #print "r2 => %s[%d:%d]" % (r2['chr'],r2['start'],r2['end'])
                    overlap.append(r2)
                    t_ov = True
            if not t_ov:
                a_ov = False
                break
        if a_ov:
            ov_chr = overlap[0]['chr']
            ov_st = min([x['start'] for x in overlap])
            ov_end = max([x['end'] for x in overlap])
            overlap_peaks.append({'chr':ov_chr,'start':ov_st,'end':ov_end,'overlap':overlap})

    intersect_name = "_".join([x.name for x in macs_tables])
    print "%s: %d overlapping regions" % (intersect_name,len(overlap_peaks))
    if (print_files):
        out = open(intersect_name + "_int.txt","w")
        print "# %s: %d overlapping regions" % (intersect_name,len(overlap_peaks))
        for p in overlap_peaks:
            print >> out, "\t".join([p['chr'],p['start'],p['end'],p['overlap']])

    for combo in it.combinations(macs_tables,len(macs_tables) - 1):
        macs_overlap(print_files,*combo)

    return overlap_peaks

def make_array_from_peaks(bind_method,macs_table,exp):
    #
    # Take coords from a table of MACS peaks, extract pileup values from
    # expts using bind_method (sum = 's', max = 'm', length-normalized sum = 'ns')
    # returns a dict of arrays by expt
    #

    pk_array = np.zeros(len(macs_table))

    
    for (i,rw) in enumerate(macs_table):
        try:
            if bind_method == 's':
                pk_array[i] = sum(exp[rw['chr']][rw['start']:rw['end']])
            elif bind_method == 'm':
                pk_array[i] = max(exp[rw['chr']][rw['start']:rw['end']])
            elif bind_method == 'ns' and rw['length'] > 0:
                pk_array[i] = (sum(exp[rw['chr']][rw['start']:rw['end']])/rw['length'])
        except KeyError:
            #print "WARNING: %s not found in chrs, skipping" % (rw['chr'],)
            continue
        except ValueError:
            #print "WARNING: %s %d:%d has no max, skipping" % (rw['chr'],rw['start'],rw['end'])
            continue

    return pk_array

        
def plot_macs_enrichment(bind_method,pks,ref_enr_expt,*secondary_enr_expts):
    """
    create an enrichment plot (similar to Bradley et. al. 20XX) from previously normalized ChIP data.  Takes ref to MACS
    table, reference chipExpt object, and list of secondary chipExpt objects
    """
    rows = []
    expts = [ref_enr_expt]
    expts.extend(secondary_enr_expts)

    pk_rws = None

    if isinstance(pks, tb.Table):
        pk_rws = pks.iterrows()
    else:
        pk_rws = pks

    for rw in pk_rws:
        r = [rw]
        for exp in expts:
            try:
                if bind_method == 's':
                    r.append(sum(exp[rw['chr']][rw['start']:rw['end']]))
                elif bind_method == 'm':
                    r.append(max(exp[rw['chr']][rw['start']:rw['end']]))
                elif bind_method == 'ns':
                    r.append((sum(exp[rw['chr']][rw['start']:rw['end']])/rw['length']))
            except KeyError:
                print "WARNING: %s not found in chrs, skipping" % (rw['chr'],)
                break
            except ValueError:
                print "WARNING: %s %d:%d has no max, skipping" % (rw['chr'],rw['start'],rw['end'])
                break
        if len(r) > 1:
            rows.append(r)


    for n,x in enumerate(rows):
        try:
            print x[1]
        except IndexError:
            print x
    rows.sort(cmp=lambda x,y: cmp(y[1],x[1]))
    plot_cmd = []
    x_ax = range(len(rows))
    marker = get_plot_marker()
    #expts.reverse()
    for (n,exp) in enumerate(expts):
        y_ax = [x[n+1] for x in rows]
        plt.plot(x_ax,y_ax,marker.next(),label=exp.name,alpha=0.6)
    plt.xlabel("Rank in " + ref_enr_expt.name)
    plt.ylabel("Binding Signal (A.U.)")
    plt.grid(True)
    plt.legend()
#        if rw.chr not in chr_refs.keys():
#            grps = [ref_enr_group].extend(secondary_enr_groups)
#            for g in grps:
#                for n in g._f_walkNodes():
#                   re.search(n.name,

    
def plot_macs_enrichment_w_overlap(bind_method,ex1_macs,ex2_macs,ex1,ex2,*secondary_enr_expts):
    """
    create an enrichment plot (similar to Bradley et. al. 20XX) from previously normalized ChIP data.  Takes ref to MACS
    table, reference chipExpt object, and list of secondary chipExpt objects
    """
    rows = []
    overlaps = []
    macs_ols = macs_overlap(False,ex1_macs,ex2_macs)
    print "macs ols: %d" % (len(macs_ols))
    print [(x['chr'],x['start'],x['end']) for x in macs_ols]
    expts = [ex1,ex2]
    expts.extend(secondary_enr_expts)
    for n,rw in enumerate(ex1_macs.iterrows()):
        r = [rw]
        for o in macs_ols:
            #print [o,rw]
            print "%d o: %s rw: %s" % (n,o['chr'],rw['chr'])
            if overlap(rw,o):
                print "overlap: %s %s %s <-> %s %s %s" % (o['chr'],o['start'],o['end'],rw['chr'],rw['start'],rw['end'])
                overlaps.append(n)
        for exp in expts:
            try:
                if bind_method == 's':
                    r.append(sum(exp[rw['chr']][rw['start']:rw['end']]))
                elif bind_method == 'm':
                    r.append(max(exp[rw['chr']][rw['start']:rw['end']]))
                elif bind_method == 'ns':
                    r.append((sum(exp[rw['chr']][rw['start']:rw['end']])/rw['length']))
            except KeyError:
                print "WARNING: %s not found in chrs, skipping" % (rw['chr'],)
                break
            except ValueError:
                print "WARNING: %s %d:%d has no max, skipping" % (rw['chr'],rw['start'],rw['end'])
                break
        if len(r) > 1:
            rows.append(r)

    for n,x in enumerate(rows):
        try:
            print x[1]
        except IndexError:
            print x
    rows.sort(cmp=lambda x,y: cmp(y[1],x[1]))
    plot_cmd = []
    x_ax = range(len(rows))
    marker = get_plot_marker()
    #expts.reverse()
    for (n,exp) in enumerate(expts):
        y_ax = [x[n+1] for x in rows]
        plt.plot(x_ax,y_ax,marker.next(),label=exp.name,alpha=0.6)
    plt.vlines(overlaps,ymin=plt.gca().get_ylim()[0],ymax=plt.gca().get_ylim()[1],colors='r',linestyles='solid',alpha=0.5,label=ex1.name + "-" + ex2.name + "Peak Overlap")
    plt.xlabel("Rank in " + ex1.name)
    plt.ylabel("Binding Signal (A.U.)")
    #plt.grid(True)
    plt.legend()

def plot_genome_region(chr,start,end,style,*expts):
    _plot_genome_region(chr,start,end,style,False,*expts)

def plot_genome_region_rel(chr,start,end,style,*expts):
    _plot_genome_region(chr,start,end,style,True,*expts)


def _plot_genome_region(chr,start,end,style,rel,*expts):
    """
    Plot base-by-base binding score for a defined genome region from an arbitrary number of ChIP experiment arrays

    can inherit style for (e.g.) stacking plots in subplots 

    y_axis fixes y_axis
    """
    
    if not style:
        style = get_plot_color()

    print expts
    for exp in expts:
        chr_slice = exp[chr][start:end]
        y_ax = exp[chr][start:end]
        if rel:
            y_ax = y_ax / np.max(y_ax)
        x_ax = np.arange(start,end)
        #print y_ax.shape
        #y2 = np.zeros(y_ax.shape[0])
        #print y2.shape
        #plt.fill_between(x_ax,y_ax,np.zeros(len(y_ax)).T,style.next(),label=exp.name,alpha=0.5)
        st = style.next()
        fill_col = st[0]
        plt.plot(x_ax,y_ax,st,label=exp.name)
        plt.fill_between(x_ax,y_ax,color=fill_col,label=exp.name,alpha=0.5)
    plt.xlabel(chr)
    plt.ylabel("Binding Signal (A.U.)")
    #plt.grid(True)
    plt.legend()
    

def stacked_plots(chr,start,end,h_regions=[],gff=None,rel_plot=[],*expts):
    style=get_line_style()
    ymax = max([(max(ex[chr][start:end])*1.2) for ex in expts])
    ylimits = [(0,ymax)] * len(expts)
    for rp in rel_plot:
        if isinstance(rp,tuple):
            for j in rp:
                ylimits[j] = (0,max([(max(expts[i][chr][start:end])*1.2) for i in rp]))
        else:
            ylimits[rp] = (0,(max(expts[rp][chr][start:end])*1.2))
    print ylimits

    if (gff):
        genes = gff.get_region_annot(chr,start,end,"gene")['gene']
    
    print ymax
    
    for (i,x) in enumerate(expts):
        #shape = (len(expts),1,i+1)#int(str(len(expts)) + "1" + str(i + 1))
        plt.subplot(len(expts),1,i+1)
        
        plotted_genes = []
        plot_genome_region(chr,start,end,style,x)
        plt.ylim(ylimits[i])
        for h in h_regions:
            plt.axvspan(h['start'], h['end'], facecolor='k', alpha=0.2)

        if (gff):
            for g in genes:
                pad = 0
                g['chr'] = '-'
                for gn in plotted_genes:
                    if overlap(g,gn) != (0,0):
                        pad = gn['pad'] + .11
                plt.axvspan(g['start'], g['end'], ymin=0.025+pad,ymax=0.125+pad,facecolor='g', alpha=0.2)
                annot_coords = ((g['start'] + g['end'])/2.0,(0.075+pad)*ymax)
                gene_name_direction = g['name'] + " >>>"
                if g['strand'] == '-':
                    gene_name_direction = '<<< ' + g['name']
                plt.annotate(gene_name_direction,xy=annot_coords,xytext=(0,20),textcoords='offset points',
                             arrowprops=None, ha='center', va='center')
                    
                g['pad'] = pad
                plotted_genes.append(g)
                print plotted_genes

def pks_to_dict(pks,chr_mod=""):
    new = []
    for x in pks:
        n_pk = {'chr':chr_mod + x['chr'],
                'start':x['start'],
                'end':x['end'],
                'summit':x['summit'],
                'length':(x['end'] - x['start'])}
        new.append(n_pk)
    return new

def get_pks_coords(sp1,sp2,pks,aln,verbose=0):
    alns = []
    no_alns = []
    split_alns = []
    aligned_sp1_pks = []
    if verbose >= 1:
        print "getting peaks for %s -> %s..." % (sp1,sp2)
    for (i,x) in enumerate(pks):
        try:
            if ((verbose >= 1) and (i % 500 == 0)):
                print "%d peaks processed" % (i,)
            pk_alns = aln.getSliceSpeciesCoord(sp1,x['chr'],x['start'],x['end'])
            if x['end'] - x['start'] == 0:
                no_alns.append(x)
            else:
                if (len(pk_alns) > 1):
                    split_alns.append((x,pk_alns))
                    if (verbose >= 2):
                        print "WARNING: %s spans %d alignment blocks; splitting" % (str(x),len(pk_alns))
                for pk_aln in pk_alns:
                    sp2_pk = {'chr':pk_aln[sp2].parent,
                              'start':pk_aln[sp2].parentStart,
                              'end':pk_aln[sp2].parentEnd,
                              'strand':pk_aln[sp2].parentStrand}
                    alns.append(sp2_pk)
                aligned_sp1_pks.append(x)
        except al.AlignError as a:
            if verbose >= 2:
                print "AlignError for %s:" % (str(x),)
                print a
            no_alns.append(x)
            continue
    if verbose >= 1:
        print "%s -> %s : %d total peaks, %d peaks aligned, %d split peaks, %d peaks not aligned" % (sp1,sp2,len(pks),len(alns),len(split_alns),len(no_alns))
    return (alns, no_alns, aligned_sp1_pks)
   

def overlap(r,r2):
    #
    # take two peak-like objects (must have 'chr', 'start', 'end' fields); return 
    # tuple of endpoints of overlap if they overlap, return (0,0) otherwise
    #
    ol_st = None
    ol_end = None

    if (r['chr'] == r2['chr']):
        if ((r['start'] <= r2['end']) and (r['start'] >= r2['start'])):
            ol_st = r['start']
        elif ((r2['start'] >= r['start']) and (r2['start'] <= r['end'])):
            ol_st = r2['start']
        
        if ((r['end'] <= r2['end']) and (r['end'] >= r2['start'])):
            ol_end = r['end']
        elif ((r2['end'] <= r['end']) and (r2['end'] > r['start'])):
            ol_end = r2['end']
        
    if ol_st and ol_end:
        return (ol_st,ol_end)
    else:
        return (0,0)

def peak_union(window,join_cut,names,filter,*pks_objects):
    """
    # 
    # window = look for overlaps in region +/- window from peak boundaries
    #
    # join_cut = fraction length of a peak needed to join; 0 = any overlap -> join
    #
    """
    start_union = []
    fin_union = []
    fin_intersect = []
    if len(names) == None:
        try:
            names = [x.name for x in pks_objects]
        except AttributeError:
            names = [""]*len(pks_objects)
    tot_over = 0
    all_overlaps = []
    expt_overlap = {} # split overlapping peaks up by experiment
    window = int(window)
    if (window >= 0):
        print "Peak window (peak width +/- %d)" % (window,)
    print "Starting Sets:"
    for (i,p) in enumerate(pks_objects):
        for r in p:
            pk = {'chr':r['chr'],
                  'start':(r['start'] - window),
                  'end':(r['end'] + window),
                  'length':(r['end'] + window) - (r['start'] - window),
                  'sets':[names[i]]}
            if 'summit' in r.keys():
                pk['summit'] = r['summit']
            if 'name' in r.keys():
                pk['name'] = r['name']
            start_union.append(pk)
        print "\t%s : %8d peaks" % (names[i],len(p))
    print "%d total peaks in starting set" % (len(start_union),)
    
    start_union.sort(key=(lambda x: (x['chr'],x['start'],x['end'])))
    
    while (len(start_union) > 0):
        #print "outer"
        r1 = start_union.pop(0)
        #print (i,r1)
        over=False
        ex=False
        while (len(start_union) > 0) and (not ex):
            #print "inner"
            r2 = start_union[0]
            #print r1
            #print r2
            ol = overlap(r1,r2)
            #print "overlap: ",
            #print ol
            if ((((ol != (0,0)) and (r1['length'] > 0)) and r2['length'] > 0) and 
                ((((float(ol[1]-ol[0])/r2['length']) > join_cut) or
                  ((float(ol[1]-ol[0])/r1['length']) > join_cut)))):
                #print "joined"
                
                r2 = start_union.pop(0)
                #print r2
                # keep an unaltered copy of each peak in all_overlaps
                if ((r1['sets'] != r2['sets']) and (len(r1['sets']) == 1)):
                    all_overlaps.append(cp.copy(r2))
                    if r1 not in all_overlaps:
                        all_overlaps.append(cp.copy(r1))

                r1['start'] = min(r1['start'],r2['start'])
                r1['end'] = max(r1['end'],r2['end'])
                r1['length'] = r1['end'] - r1['start']
                r1['sets'] = r1['sets'] + r2['sets']
                #r2['sets'] = r1['sets'] + r2['sets']
                tot_over += 1
                over = True
            elif (ol[1] - ol[0]) > 0:
                #print "split"
                #print ol
                r1['end'] = r2['start'] - 1
            else:
                #print "appended:"
                #print r1
                if (over):
                    tot_over += 1
                    fin_intersect.append(r1)
                fin_union.append(r1)
                ex = True
    if filter:
        #print "final isct:"
        #print fin_intersect
        #print "+++++++++++"
        flt_intersect = [x for x in fin_intersect if (len(set(x['sets'])) > 1)]
        #flt_all_overlaps = [x for x in all_overlaps if (len(set(x['sets'])) > 1)]
        for n in names:
            expt_overlap[n] = [x for x in all_overlaps if n in x['sets']]
    else:
        flt_intersect = fin_intersect
        flt_all_overlaps = all_overlaps
    print "%d in union, %d in intersect" % (len(fin_union),len(fin_intersect))
    if filter:
        print "final isct: %d before filtering, %d after" % (len(fin_intersect),len(flt_intersect))
        print "total overlapping peaks (all sets): %d" % (len(all_overlaps))#,len(flt_all_overlaps))
        print "By-Experiment Overlap Counts:"
        for (n,a) in expt_overlap.items():
            print "\t%s : %8d" % (n,len(a))
    return (fin_union,flt_intersect,expt_overlap)
    
def count_chr_dist(pks):
    chr_cnt = {}
    chr_prop = {}
    total = 0
    for c in [x['chr'] for x in pks]:
        if (c in chr_cnt.keys()):
            total += 1
            chr_cnt[c] += 1
        else:
            total += 1
            chr_cnt[c] = 1
    chr_prop = dict([(sp,float(c)/total) for (sp,c) in chr_cnt.items()])
    print "====== counts:"
    for i in sorted(chr_cnt.items()):
        print "%s => %d" % i
    print "====== proportions:"
    for i in sorted(chr_prop.items()):
        print "%s => %f" % i
    return chr_cnt,chr_prop

def map_fn_to_aln_list(coord_list,ref_sp,aln,fn,out_it,pad=0,save_alns=None,verbose=False,fill_aln=False):
    out = cp.copy(out_it)
    for (i,c) in enumerate(coord_list):
        if i % 100 == 0:
            print "%d coordinates processed (%4f%%..." % (i,(float(i)/len(out)) * 100)
        try:
            bs_aln = aln.getSliceSpeciesCoord(ref_sp,c['chr'],c['start']-pad,c['end']+pad,fill_missing=fill_aln)
            try:
                if (('strand' in c.dtype.names) and c['strand'] == "-"):
                    bs_aln[0] = bs_aln[0].revcomAln()
                elif (('strand' in c.dtype.names) and c['strand'] == "+"):
                    bs_aln[0] = bs_aln[0].revcomAln()
            except:
                if ((('strand' in c.keys()) and c['strand'] == "-")):
                    bs_aln[0] = bs_aln[0].revcomAln()
        except al.AlignError as al_err:
            if verbose:
                print "Align Error for coordinate: ",
                print c
                print al_err
            continue
        except AttributeError as att_e:
            print att_e
        #print bs_aln[0]
        if len(bs_aln) > 1:
            print "chimeric aln - skipped",
            print c
            continue
        out[i] = fn(bs_aln[0])
    return out

def patser_macs_overlap(peaks,phits):
    count = 0
    bnd_coords = []
    for (i,p) in enumerate(peaks):
        if i % 100 == 0:
            print (i,p)
            wl = phits.getWhereList('(chr == \'%s\') & ((start > %d) & (end < %d))' % (p['chr'],p['start'],p['end']))
            print len(wl)
            print "%d pks %4f%% done, %d patser hits" % (i,(float(i)/len(peaks)) * 100,len(bnd_coords))
        h = phits.getWhereList('(chr == \'%s\') & ((start > %d) & (end < %d))' % (p['chr'],p['start'],p['end']))
        count += len(h)
        bnd_coords.extend(h)
    bnd = phits.readCoordinates(bnd_coords)
    bnd['overlap'] = True
    print "updating table..."
    phits.modifyRows(bnd_coords,bnd)
    #bnd = phits.readWhere('overlap')
    print "%d peaks, %d hits, %d hits in bound regions, %d outside" % (len(peaks),len(phits),len(bnd),(len(phits) - len(bnd)))
    #pat_hit_dtype = np.dtype([('chr','S25'),('start',np.uint32),('end',np.uint32),('name','S20'),('score',np.float32),('strand','S1'),('pval',np.float32),('peak_over',np.uint8)])
    #bnd = np.array(bnd,dtype=pat_hit_dtype)
    phits.flush()
    return bnd

def get_top_n_peaks(n,peaks,tag_h5,input_h5=None):
    n_peaks = []
    pk_array = make_array_from_peaks("ns",peaks,tag_h5)
    if input_h5:
        input_array = make_array_from_peaks("ns",peaks,input_h5)
        fin_array = pk_array - input_array
        for (pe,pi) in zip(peaks,fin_array):
            n_peaks.append((pi,pe))
    else:
        for (pe,pi) in zip(peaks,pk_array):
            n_peaks.append((pi,pe))
    srt_pks = sorted(n_peaks, cmp=(lambda x,y: cmp(y[0],x[0])))
    return srt_pks[:n-1]

def fit_line_to_pts(x,y,plot=True,residuals=False):
    # fit a line y=mx+c to equal-length vectors x,y
    A = np.vstack([x, np.ones(len(x))]).T
    regr = np.linalg.lstsq(A, y)
    (m,c) = regr[0]
    if (plot):
        plt.plot(x, m*x + c, 'r', label='Fitted line')
    return (m,c)

def scatter_w_reg_line(x,y,pnt_lab,pt_col,xlab,ylab,title,ln=True,xyline=True):
    plt.plot(x,y,"."+pt_col,alpha=0.4,label=pnt_lab)
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    if ln:
        plt.yscale('log')
        plt.xscale('log')
    if xyline:
        plt.plot([plt.xlim()[0],plt.xlim()[1]],[plt.xlim()[0],plt.xlim()[1]],":k",alpha=0.5,label="x == y")
    (pm,pc) = fit_line_to_pts(x,y,plot=False)
    xpts = np.arange(plt.xlim()[0],plt.xlim()[1],abs(plt.xlim()[0] - plt.xlim()[1])/100)
    plt.plot(xpts,pm*xpts + pc,"--"+pt_col,label="Lst. Sq. (m=%f,c=%f)" % (pm,pc))
    plt.title(title)
    plt.legend(loc=2)
    return (pm,pc)

def fit_line_to_pts_w_residuals(x,y,plot=True):
    # fit a line y=mx+c to equal-length vectors x,y
    A = np.vstack([x, np.ones(len(x))]).T
    regr = np.linalg.lstsq(A, y)
    (m,c) = regr[0]
    res = regr[1]
    if (plot):
        plt.subplot(211)
        plt.plot(x, m*x + c, 'r', label='Fitted line')
        plt.subplot(212)
        plt.hist(res,bins=100)
    return (m,c,res)

def cluster_by_target(pks,gns,name_idx,sp):
    cls = {}
    for (p,g) in zip(pks,gns):
        if g[name_idx] in cls.keys():
            if sp in cls[g[name_idx]].keys():
                cls[g[name_idx]][sp].append(p)
            else:
                cls[g[name_idx]][sp] = [p]
        else:
            cls[g[name_idx]] = {sp:[p],
                                sp+"_gene" : g}
    return cls

def clusters_overlap(dm_cl,dp_cl):
    ol_cls = {}
    for g in dm_cl.keys():
        if g in dp_cl.keys():
            ol_cls[g] = cp.copy(dp_cl[g])
            ol_cls[g].update(cp.copy(dm_cl[g]))
    for g in ol_cls.keys():
        print "%s : %s dmel peaks, %s dpse peaks" % (g,len(ol_cls[g]['dmel']),len(ol_cls[g]['dpse']))
    return ol_cls

def cluster_peaks_by_feature(pks,gff,feature_type='gene',name_idx=3,sp='dmel'):
    gns = [gff.get_nearest_annot(x,5,feature_type)[feature_type][0][1] for x in pks]
    cl = cluster_by_target(pks,gns,name_idx,sp)
    return cl

def ol_genes(dm_pks,dp_pks,dm_gff,dp_gff):
    dm_gns = [dm_gff.get_nearest_annot(x,5,'gene')['gene'][0][1] for x in dm_pks]
    dp_gns = [dp_gff.get_nearest_annot(x,5,'orthologous_to')['orthologous_to'][0][1] for x in dp_pks]
    
    dm_cl = cluster_by_target(dm_pks,dm_gns,3,'dmel')
    dp_cl = cluster_by_target(dp_pks,dp_gns,1,'dpse')
    
    ol_cl = clusters_overlap(dm_cl,dp_cl)
    
    print "\n%d dmel clusters, %d dpse clusters, %d gene overlaps" % (len(dm_cl.keys()),len(dp_cl.keys()),len(ol_cl.keys()))
    
    vn.venn2(subsets=((len(dm_cl.keys()) - len(ol_cl.keys())),(len(dp_cl.keys()) - len(ol_cl.keys())),len(ol_cl.keys())))
    
    return dm_cl,dp_cl,ol_cl

def plot_by_diff_rank(plot,diff_m,revxy,peaks,expt_x,expt_y,*extra):
    #
    # make an enrichment plot ranked by (expt_x - expt_y)
    #
    # plot toggles plotting on / off
    #
    # diff_m="rel" means relative difference, diff_m="abs" means absolute (*NOT* absolute value) diff
    #
    # revxy sets all expt y values negative, so they appear below the x axis (True|False)
    rows = []
    expts = [expt_x,expt_y]
    expts.extend(extra)
    dm = ""
    
    row_itr = None
    
    for rw in peaks:
        
        r = None
        
        try:
            x_nsum = sum(expt_x[rw['chr']][rw['start']:rw['end']])/rw['length'] 
            y_nsum = sum(expt_y[rw['chr']][rw['start']:rw['end']])/rw['length']
            diff = x_nsum - y_nsum
            
            if (x_nsum == 0) and (y_nsum == 0):
                continue
            else:
                r = [rw]
            
            if diff_m == 'rel':
                dm = "Relative"
                diff /= sum(expt_x[rw['chr']][rw['start']:rw['end']])/rw['length']
            elif diff_m == 'abs':
                dm = "Absolute"
            elif diff_m == 'abs_v':
                diff = abs(diff)


            r.append(diff)
                
        except KeyError:
            print "WARNING: %s not found in chrs, skipping" % (rw['chr'],)
            print rw
            break
        except ValueError:
            print "WARNING: %s %d:%d has no max %s, skipping" % (rw['chr'],rw['start'],rw['end'],"diff")
            break

        for exp in expts:
            try:
                r.append((sum(exp[rw['chr']][rw['start']:rw['end']])/rw['length']))
            except KeyError:
                print "WARNING: %s not found in chrs, skipping" % (rw['chr'],)
                print rw
                break
            except ValueError:
                print "WARNING: %s %d:%d has no max, skipping" % (rw['chr'],rw['start'],rw['end'])
                break
        if len(r) > 1:
            rows.append(r)
    
    rows.sort(cmp=lambda x,y: cmp(y[1],x[1]))
    plot_cmd = []
    x_ax = range(len(rows))
    marker = get_plot_marker()
    #expts.reverse()
    y_ax = [x[1] for x in rows]
    if (plot):
        plt.plot(x_ax,y_ax,marker.next(),label=dm + " Difference",alpha=0.6)
    for (n,exp) in enumerate(expts):
        mlt = 1
        if revxy and (n==1):
            mlt = -1
        y_ax = [mlt*(x[n+2]) for x in rows]
        if plot:
            plt.plot(x_ax,y_ax,marker.next(),label=exp.name,alpha=0.6)
    if plot:
        plt.xlabel("Rank in %s - %s (%s Difference)" % (expt_x.name,expt_y.name,dm))
        plt.ylabel("Binding Signal (A.U.)")
        plt.grid(True)
        plt.legend()
    
    return rows

def diff_enr_plots(color_map,mtx_max,peaks,mark_ranks,x,y,*others):
    sorted_rows = plot_by_diff_rank(False,'abs',True,peaks,x,y,*others)
    mtx = []
    order = [0] + range(2,len(others) + 2) + [1] # note that x is always on top, y always on bottom
    for i in order:
        diff_srt = [rw[i+2] for rw in sorted_rows]
        mtx.append(diff_srt)
    diff_mtx = np.array(mtx)
    plt.imshow(diff_mtx,cmap=color_map,aspect='auto',interpolation='none',vmax=mtx_max)
    if len(mark_ranks) > 0:
        plt.vlines(mark_ranks,plt.ylim()[0],plt.ylim()[1],color="white",linewidth=2,linestyle="dashed")

def blast_seq(query,genome=DMEL_GENOME,name="query_seq",top=5):
        all_hits = []
        fa_seq = ">%s\n%s\n" % (name,query)
        temp_q = tmp.NamedTemporaryFile(suffix=".fa")
        temp_q_name = temp_q.name
        print >> temp_q, fa_seq
        temp_q.flush()
        blast_proc = sub.Popen([BLAST_EXE,genome,temp_q_name,"-mformat","3","-sort_by_highscore","-sort_by_subjectlength"],stdout=sub.PIPE,stderr=sub.PIPE)
        (out,err) = blast_proc.communicate()
        temp_q.close()
        #print out
        #print err
        hits = sorted([x.strip().split("\t") for x in out.split("\n") if (len(x) > 1) and (x[0] != "#")],cmp=(lambda z,y: cmp(float(y[4]),float(z[4]))))
        #print hits
        comments = [x for x in out.split("\n") if (len(x) > 1) and (x[0] == "#")]
        for bl_ln_sp in hits[:top]:
            bl = {'chr':bl_ln_sp[1],
                  'start':int(bl_ln_sp[20]),
                  'end':int(bl_ln_sp[21]),
                  'length':int(bl_ln_sp[21]) - int(bl_ln_sp[20]),
                  'strand':"+",
                  'name':name}
            if bl_ln_sp[16] == "-1":
                bl['strand'] = "-"
            all_hits.append(bl)
        return all_hits

############
#
# Description Classes for pytables
#

class MacsPeak(tb.IsDescription):

    chr = tb.StringCol(32)
    start = tb.UInt32Col()
    end = tb.UInt32Col()
    length = tb.UInt16Col()
    summit = tb.UInt32Col()
    tags = tb.UInt32Col()
    log_pval = tb.Float32Col()
    fold_enrichment = tb.Float32Col()
    log_qval = tb.Float32Col()

class gff_record(tb.IsDescription):

    start = tb.UInt32Col()
    end = tb.UInt32Col()
    name = tb.StringCol(10)
    strand = tb.StringCol(1)
    FBid = tb.StringCol(11)
    Parent = tb.StringCol(11)
    
############
#
# Other Classes
#

class chipExpt():
    def __init__(self,group,name="NoName",chrs=DMEL_CHRS,chr_ext='ext'):
        self.chr_dict = {}
        self.group = group
        self.name = name
        arrays = group._v_children.items()
        if chrs:
            for ch in chrs:
                for (name,node) in arrays:
                    if re.search(chr_ext + "\w*_(chr)?"+ch+"$",name):
                        #print (ch,name,node.name)
                        self.chr_dict[ch] = node
                        break
        else:
            for n in group._f_walkNodes():
                ch_s = re.search("chr(.+)",n.name)
                self.chr_dict[ch_s.group(1)] = n
                break
                    
    def __getitem__(self,item):
        return self.chr_dict[item]

    def write_wig(self,chr_convert=None,span=10,outfile=None):
        if outfile == None:
            outfile = self.name + ".wig"
        wigout = open(outfile,"w")
        print >> wigout, "track type=wiggle_0 name=%s description=%s visibility=full autoScale=off maxHeightPixels=100:50:20" % (self.name,self.name)
        for (chrom,chr_arr) in self.chr_dict.items():
            print "Writing %s..." % (chrom,),
            name = chrom
            if chr_convert != None:
                name = chr_convert[chrom]
            print >> wigout, "variableStep chrom=%s span=%d" % (name,span)
            for i in np.arange(1,len(chr_arr),step=span):
                print >> wigout, "%d\t%f" % (i,np.sum(chr_arr[i:i+span])/float(span))
            print "done!"
        wigout.close()
                
class chipExpt_wig():
    def __init__(self,file,name="NoName"):#,chrs=DMEL_CHRS,chr_ext='ext'):
        self.chr_dict = {}
        self.file = file
        self.name = name
        self._parse_file()

    def _parse_file(self):
        f = open(self.file)
        curseq = None
        start = None
        step = None
        fixed = False
        prev_coord = None
        vals = []
        for ln in f:
            if re.match('track',ln):
                if self.name == "NoName":
                    name_match = re.match('name=(\S+)',ln)
                    self.name = name_match.group(1)
            elif re.match('variableStep',ln):
                if curseq:
                    self.chr_dict[curseq] = self._wigArray(vals,start,step)
                    vals=[]
                    start = None
                    step = None
                vs_match = re.search('chrom=.*?chr(\S+).+span=(\S+)',ln)
                curseq = vs_match.group(1)
                step = int(vs_match.group(2))
            elif re.match('fixedStep',ln):
                fixed = True
                if curseq:
                    self.chr_dict[curseq] = self._wigArray(vals,start,step)
                    vals=[]
                    start = None
                    step = None
                vs_match = re.search('chrom=(chr)+(\S+).+start=(\S+).+step=(\S+)',ln)
                curseq = vs_match.group(2)
                start = int(vs_match.group(3))
                step = int(vs_match.group(4))

            elif fixed:
                ln_sp = ln.strip()
                #print ln_sp
                #if not start:
                #    start = int(float(ln_sp[0]))
                #if (prev_coord) and (int(ln_sp[0]) - prev_coord != step):
                #    print "WARNING: Jumped coords %s: %d -> %d" % (curseq,prev_coord,ln_sp_0)
                vals.append(np.float64(ln_sp))
            else:
                ln_sp = ln.strip().split()
                if not start:
                    start = int(float(ln_sp[0]))
                if (prev_coord) and (int(ln_sp[0]) - prev_coord != step):
                    print "WARNING: Jumped coords %s: %d -> %d" % (curseq,prev_coord,ln_sp_0)
                vals.append(np.float64(ln_sp[1]))
        if curseq:
            self.chr_dict[curseq] = self._wigArray(vals,start,step)

                
                
    def __getitem__(self,item):
        return self.chr_dict[item]

    class _wigArray:

        def __init__(self,ls,start,step,smooth_win=5,name="NoName"):
            self.ar = np.array(ls)
            self.name = name
            self.start = start
            self.step = step
            self.smoothwin = step*smooth_win

        def __getslice__(self,start,end,smoothed=True):
            orig_start = start
            orig_end = end
            if end - start > self.smoothwin:
                end = end + int(self.smoothwin/2.0)
                start = start - int(self.smoothwin/2.0)
            start_idx = int((start - self.start)/self.step) + 1
            start_remainder = ((start_idx*self.step) + self.start) - start
            end_idx = int((end-self.start)/self.step)
            end_remainder = end - ((end_idx*self.step) + self.start)
            compressed_vals = self.ar[start_idx:end_idx]
            fin_array = []
            for v in compressed_vals:
                fin_array.extend([float(v)]*self.step)
            fin_array = [compressed_vals[0]]*start_remainder + fin_array + [compressed_vals[-1]] * end_remainder
            if end - start > self.smoothwin:
                fin_array = [sum(fin_array[n:n+self.smoothwin])/(self.smoothwin) for n in range(0,orig_end - orig_start)]
            #coords = np.arange(start_idx*self.step + self.start,end_idx*self.step + self.start,step=self.step)
            return np.array(fin_array)

class GFFtable():
    def __init__(self,gff_file,h5_file,chrs,*entry_types):
        self.entry_types = entry_types
        self.file = gff_file
        self.chrs = chrs
        self.h5 = None
        if (not os.path.exists(h5_file)):
            self.h5 = tb.openFile(h5_file,"w")
            self._parse_gff(gff_file,entry_types,chrs)
        else:
            self.h5 = tb.openFile(h5_file,"r+")

    def __del__(self):
        self.h5.close()

    def get_region_annot(self,chr,start,end,*entry_types):
        """
        get annotations for a specified region, return dict w/ keys as entry type, pointing to list of entries
        """
        e_dict = {}
        print entry_types
        for e in entry_types:
            e_dict[e] = []
            tab = self.h5.getNode("/" + e + "/" + "chr" + chr)
            query = "(((start >= %d) & (start <= %d)) | ((end >= %d) & (end <= %d)))" % (start,end,start,end)
            for x in tab.where(query):
                g = {'name':x['name'],
                     'start':x['start'],
                     'end':x['end'],
                     'strand':x['strand'],
                     'FBid':x['FBid'],
                     'Parent':x['Parent']}
                print x['name']
                print x['start']
                e_dict[e].append(g)
            print [(x['name'],x['start']) for x in e_dict[e]] 
        return e_dict

    def get_nearest_annot(self,peak,n_nearest,*entry_types):
        pos = int((peak['start'] + peak['end'])/2) # hack for now - just take exact middle of region - could fix to get summit values later (2/8/13)
        e_dict = {}
        for e in entry_types:
            e_dict[e] = []
            tab = self.h5.getNode("/" + e + "/" + "chr" + peak['chr'])
            if len(tab) == 0:
                break
            dist_ar = np.zeros(len(tab),dtype=np.int32)
            dist_ar[:] = map((lambda x: self.start_dist(x,pos)),tab)
            #print peak
            #print dist_ar
            d_max = np.max(abs(dist_ar))
            idx_min = np.argmin(abs(dist_ar))
            e_dict[e] = [(dist_ar[idx_min],tab[idx_min])]
            dist_ar[idx_min] = d_max
            for i in range(0,n_nearest - 1):
                idx_min = np.argmin(abs(dist_ar))
                e_dict[e].append((dist_ar[idx_min],tab[idx_min]))
                dist_ar[idx_min] = d_max
        return e_dict

    def start_dist(self,entry,posn):
        if entry['strand'] == '+':
            return posn - entry['start']
        else:
            return posn - entry['end']

    def gene_location(self,gene_name,gene_list=False,entry_type="gene"):
        hits = []
        for chr_node in self.h5.iterNodes("/" + entry_type):
            chr_name_re = re.search("chr(.+)",chr_node.name)
            chr_name = chr_name_re.group(1)
            hiterator = chr_node.where("name == \"%s\"" % (gene_name,))
            for h in hiterator:
                loc={'start' : h['start'],
                     'end' : h['end'],
                     'strand' : h['strand'],
                     'chr' : chr_name,
                     'name' : gene_name,
                     'FBid' : h['FBid'],
                     'Parent' : h['Parent']}
                hits.append(loc)
        if gene_list:
            return hits

        if len(hits) == 0:
            print "WARNING: no hits for gene name %s!" % (gene_name,)
            return []
        elif len(hits) > 1:
            print "WARNING: multiple hits for gene name %s!" % (gene_name,)

        return hits[0]

    def get_all_entry_type(self,entry_type):
        hits = []
        for chr_node in self.h5.iterNodes("/" + entry_type):
            chr_name_re = re.search("chr(.+)",chr_node.name)
            chr_name = chr_name_re.group(1)
            if chr_name not in self.chrs:
                continue
            for h in chr_node.iterrows():
                loc={'start' : h['start'],
                     'end' : h['end'],
                     'strand' : h['strand'],
                     'chr' : chr_name,
                     'name' : h['name'],
                     'FBid' : h['FBid'],
                     'Parent' : h['Parent']}
                hits.append(loc)
        return hits


    def _parse_gff(self,gff_file,entry_types,chrs):

        t_dict = {}
        for e in entry_types:
            e_grp = self.h5.createGroup("/",e,e + " annotations")
            t_dict[e] = {}
            for c in chrs:
                t_dict[e]["chr"+c] = self.h5.createTable(e_grp,"chr"+c,gff_record,e + " annotations, " + c)
                g = t_dict[e]["chr"+c]
                g.cols.start.createIndex()
                g.cols.end.createIndex()

        for l in open(gff_file):
            if re.match("#",l):
                continue
            m = re.search("^(\S+)\s+\S+\s+(\S+)\s+(\S+)\s+(\S+)\s+\S+\s+(\S)\s+\S+\s+(.+)",l)
            if not m:
                print "WARNING: no line match for:\n%s" % (l,)
                continue
            c = m.group(1) 
            e_type = m.group(2)
            if (not e_type in entry_types) or (not c in chrs):
                continue
            entry = t_dict[e_type]["chr"+c].row
            entry['start'] = m.group(3)
            entry['end'] = m.group(4)
            entry['strand'] = m.group(5)
            info_string = m.group(6)
            #print info_string
            id_m = re.search("ID=([^;]+);",info_string)
            name_m = re.search("Name=([^;]+);",info_string)
            parent_m = re.search("Parent=([^;]+);",info_string)
            to_name_m = re.search("to_name=([^;]+);",info_string)
            gene_id_m = re.search("gene_id \"([^;]+)\";",info_string)
            #transcript_id_m = re.search("transcript_id \"([^;]+)\";",info_string)

            if id_m:
                entry['FBid'] = id_m.group(1)
            if name_m:
                entry['name'] = name_m.group(1)
            if parent_m:
                entry['Parent'] = parent_m.group(1)
            if to_name_m:
                entry['Parent'] = to_name_m.group(1)
            if gene_id_m:
                entry['name'] = gene_id_m.group(1)
            #if transcript_id_m:
            #    entry['transcript_id'] = transcript_id_m.group(1)
            #print "%s %s %s %s %s %s %s" % (c,entry['start'],entry['end'],entry['strand'],entry['FBid'],entry['name'],entry['Parent'])
            entry.append()

        #print e    
        for e in entry_types:
            for c in chrs:
                t_dict[e]["chr"+c].flush()
            
