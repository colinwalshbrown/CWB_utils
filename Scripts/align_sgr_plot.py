#!/usr/bin/env python
"""
Get a slice from whole genome alignment, make plot using sgr for each sp
"""

import sys
import SGR
import alignment
import matplotlib.pyplot as plt


def main(args):

    if len(args) < 9:
        print "usage: align_sgr_plot.py <sp1_name> <sp1_sgr> <sp2_name> <sp2_sgr> <align> <chr> <start> <end> <plot order 0:left=top,1:right=top>"
        sys.exit(1)

    plotSGRslice(args[0],args[1],args[2],args[3],args[4],args[5],int(args[6]),int(args[7]),int(args[8]))
        
    
def plotSGRslice(sp1_name,sp1_sgr,sp2_name,sp2_sgr,align,chr,start,end,order):
    sgr1 = SGR.SGR(sp1_sgr,reindex=False)
    sgr2 = SGR.SGR(sp2_sgr,reindex=False)
    aln = alignment.WholeGenomeAlign(align,format="FSA")
    aln_slice = aln.getSliceSpeciesCoord(sp1_name,chr,start,end)[0]
    if aln_slice.alignSeqs[sp1_name].parentStrand == "-":
        aln_slice = aln_slice.revcomAln()
    sgr1_slice = sgr1.SGRslice(chr,start,end)
    sp2_chr = aln_slice.alignSeqs[sp2_name].parent
    sp2_start = aln_slice.alignSeqs[sp2_name].parentStart
    sp2_end = aln_slice.alignSeqs[sp2_name].parentEnd
    sp2_strand = aln_slice.alignSeqs[sp2_name].parentStrand
    print (sp2_chr,sp2_start,sp2_end,sp2_strand)
    sgr2_slice = sgr2.SGRslice(sp2_chr,sp2_start,sp2_end)

    print aln_slice.length

    counts = {sp1_name : [],
              sp2_name : []}

    sp1_alnseq = aln_slice.alignSeqs[sp1_name].seq
    sp2_alnseq = aln_slice.alignSeqs[sp2_name].seq

    sp1_sgr_idx = 0
    sp2_sgr_idx = 0

    sp1_vals = expand_sgr(sgr1_slice['vals'],start,end,"+",sgr1_slice['step'])
    sp2_vals = expand_sgr(sgr2_slice['vals'],sp2_start,sp2_end,sp2_strand,sgr2_slice['step'])

    sp1_cur = sp1_vals[sp1_sgr_idx]
    sp2_cur = sp2_vals[sp2_sgr_idx]

    for x in range(0,aln_slice.length):
        print "%s 1:%s ix:%s cur:%s 2:%s ix:%s cur:%s" % (x,sp1_alnseq[x],sp1_sgr_idx,sp1_cur,sp2_alnseq[x],sp2_sgr_idx,sp2_cur) 
        
        if sp1_alnseq[x] == "-":
            counts[sp1_name].append(sp1_cur)
        else:
            if sp1_sgr_idx < len(sp1_vals) -1:
                sp1_sgr_idx += 1
            sp1_cur = sp1_vals[sp1_sgr_idx]
            counts[sp1_name].append(sp1_cur)
            
        if sp2_alnseq[x] == "-":
            counts[sp2_name].append(sp2_cur)
        else:
            if sp2_sgr_idx < len(sp2_vals) - 1:
                sp2_sgr_idx += 1
            sp2_cur = sp2_vals[sp2_sgr_idx]
            counts[sp2_name].append(sp2_cur)

        #print counts

    """
        sp1_unaln = aln_slice.alignSeqs[sp1_name].getUnalignedCoord(x)
        sp1_pcoor = start + sp1_unaln
        if sp1_pcoor % sgr1_slice['step'] == 0:
            sp1_sgr_idx += 1
        sp2_unaln = aln_slice.alignSeqs[sp2_name].getUnalignedCoord(x)
        if sp2_strand == "-":
            sp2_pcoor = sp2_end - sp2_unaln
            if sp2_pcoor % sgr2_slice['step'] == 0:
                sp2_sgr_idx -= 1
        else:
            sp2_pcoor = sp2_start + sp2_unaln
            if sp2_pcoor % sgr2_slice['step'] == 0:
                sp2_sgr_idx += 1

        if (sp2_sgr_idx <= 0):
            sp2_sgr_idx = 0
        elif (sp2_sgr_idx >= len(sgr2_slice['vals'])):
            sp2_sgr_idx = len(sgr2_slice['vals']) - 1
        
        if (sp1_sgr_idx >= len(sgr1_slice['vals'])):
            sp1_sgr_idx = len(sgr1_slice['vals']) - 1         
         
        counts[sp1_name].append(sgr1_slice['vals'][sp1_sgr_idx])
        counts[sp2_name].append(sgr2_slice['vals'][sp2_sgr_idx])
    """
    print order
    if order == 0:
        plt.figure(1,figsize=(5,6))
        plt.subplot(211)
        if (max(counts[sp1_name])<20):
            plt.ylim(0,20)
        plt.plot(counts[sp1_name],lw=3)
        plt.ylabel("Normalized Tags",fontsize=14)
        plt.xlabel("%s:%s-%s" % (chr,start,end))
    #plt.text(1,1,sp1_name,va="top",ha="right",size=14,transform=ax.transAxes)
        plt.grid(True)

        plt.subplot(212)
        plt.plot(counts[sp2_name],'r',lw=3)
        if (max(counts[sp2_name])<20):
            plt.ylim(0,20)
        plt.ylabel("Normalized Tags",fontsize=14)
        plt.xlabel("%s:%s-%s" % (sp2_chr,sp2_start,sp2_end))
    #plt.text(1,1,sp2_name,va="top",ha="right",size=14,transform=ax.transAxes)
        plt.grid(True)
    else:
        plt.figure(1,figsize=(5,6))
        plt.subplot(211)
        if (max(counts[sp2_name])<20):
            plt.ylim(0,20)
        plt.plot(counts[sp2_name],lw=3)
        plt.ylabel("Normalized Tags",fontsize=14)
        plt.xlabel("%s:%s-%s" % (sp2_chr,sp2_start,sp2_end))
    #plt.text(1,1,sp2_name,va="top",ha="right",size=14,transform=ax.transAxes)
        plt.grid(True)

        plt.subplot(212)
        plt.plot(counts[sp1_name],'r',lw=3)
        if (max(counts[sp1_name])<20):
            plt.ylim(0,20)
        plt.ylabel("Normalized Tags",fontsize=14)
        plt.xlabel("%s:%s-%s" % (chr,start,end))
    #plt.text(1,1,sp1_name,va="top",ha="right",size=14,transform=ax.transAxes)
        plt.grid(True)
        

    plotfile = "%s:%s:%s-%s-%s:%s:%s-%s.svg" % (sp1_name,chr,start,end,sp2_name,sp2_chr,sp2_start,sp2_end)
    alnfile = open("%s:%s:%s-%s-%s:%s:%s-%s.clw" % (sp1_name,chr,start,end,sp2_name,sp2_chr,sp2_start,sp2_end),"w+")
    
    aln_slice.PrintClustal(60,outHandle=alnfile)

    plt.savefig(plotfile,format='svg')
    plt.show()

def expand_sgr(values,start,end,strand,step):
    print "%s %s %s %s %s" % (len(values),start,end,strand,step)
    exp_vals = []
    idx = 0
    val = values[idx]
    for x in range(start,end):
        print "%s %s %s" % (x, idx, val)
        exp_vals.append(val)
        if x % step == 0:
            if idx < (len(values) - 1):
                idx += 1
            val = values[idx]
    print len(exp_vals)
    if strand == "-":
        exp_vals.reverse()
    return exp_vals
            

if __name__ == "__main__":
    main(sys.argv[1:])
