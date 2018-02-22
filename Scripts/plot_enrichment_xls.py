#!/usr/bin/env python
"""
Get a slice from whole genome alignment, make plot using sgr for each sp
"""

import re
import sys
import SGR
import alignment
import matplotlib.pyplot as plt
#import scipy.stats.stats as st
import numpy

def main(args):

    if len(args) < 8:
        print "usage: align_sgr_plot.py <sgr1> <sgr1_Name> <sgr2> <sgr2_name> <In_sgr> <xls> <output_base> <0=plot_rank,1=plot_scatter,2=histogram>"
        sys.exit(1)
        
    sgr1 = SGR.SGR(args[0])
    sgr2 = SGR.SGR(args[2])
    sgrIn = SGR.SGR(args[4])
    outbase = args[6]
    plot_type = args[7]

    xls_lines = []
    xls_entries = []

    for line in open(args[5]):
        if re.match("#",line):
            continue
        spl = line[:-1].split()
        entr = {'chr':spl[1],
                'start':int(spl[2]),
                'end':int(spl[3])}
        xls_lines.append(entr)

    for x in xls_lines:
        try:
            sgr1_slice = sgr1.SGRslice(x['chr'],x['start'],x['end'])
            sgr2_slice = sgr2.SGRslice(x['chr'],x['start'],x['end'])
            sgrIn_slice = sgrIn.SGRslice(x['chr'],x['start'],x['end'])
        except:
            print x
            continue
        sgr1_sum = sum(sgr1_slice['vals'])
        sgr2_sum = sum(sgr2_slice['vals'])
        sgrIn_sum = sum(sgrIn_slice['vals'])

        x['sgr1_diff'] = sgr1_sum - sgrIn_sum
        if sgrIn_sum > 0:
            x['sgr1_enr'] = sgr1_sum / sgrIn_sum
        else:
            x['sgr1_enr'] = sgr1_sum

        x['sgr2_diff'] = sgr2_sum - sgrIn_sum
        if sgrIn_sum > 0:
            x['sgr2_enr'] = sgr2_sum / sgrIn_sum
        else:
            x['sgr2_enr'] = sgr2_sum
        xls_entries.append(x)

    diff_sgr1_sorted = sorted(xls_entries,key=(lambda x:x['sgr1_diff']))
    diff_sgr2_sorted = sorted(xls_entries,key=(lambda x:x['sgr2_diff']))

    diff_sgr1_sorted.reverse()
    diff_sgr2_sorted.reverse()
    
#    print diff_sgr1_sorted
#    print diff_sgr1_sorted

    enr_sgr1_sorted = sorted(xls_entries,key=(lambda x:x['sgr1_enr']))
    enr_sgr2_sorted = sorted(xls_entries,key=(lambda x:x['sgr2_enr']))

    enr_sgr1_sorted.reverse()
    enr_sgr2_sorted.reverse()

    plt.figure(1,figsize=(6,6))

    if plot_type == '0':
        # draw and save a ranked-enrichment plot
        plt.plot([x['sgr1_diff'] for x in diff_sgr1_sorted],'rs',alpha=0.5,label=args[1])
        plt.plot([x['sgr2_diff'] for x in diff_sgr1_sorted],'bs',alpha=0.5,label=args[3])
        plt.ylabel("Normalized, Input-Subtracted Tags")
        plt.xlabel(args[1] + " Rank")
        #plt.ylim((-100,700))
        plt.legend(numpoints=1,scatterpoints=1)
        plt.grid(True)
    elif plot_type == '1':
        # draw and save a scatterplot
        plt.plot([x['sgr1_diff'] for x in enr_sgr1_sorted],[x['sgr2_diff'] for x in enr_sgr1_sorted],'oc',alpha=0.7)
        plt.ylabel(args[3] + ": Normalized, Input-Subtracted Tags")
        plt.xlabel(args[1] + ": Normalized, Input-Subtracted Tags")
        plt.grid(True)
    elif plot_type == '2':
        # draw and save a histogram
        plt.hist([x['sgr1_diff'] for x in enr_sgr1_sorted],bins=50,alpha=0.5,label=args[1])
        plt.hist([x['sgr2_diff'] for x in enr_sgr1_sorted],bins=25,alpha=0.5,label=args[3])
        plt.legend(numpoints=1,scatterpoints=1)
        diff1_out = open(outbase + "_sgr1_diff_list.txt","w")
        diff2_out = open(outbase + "_sgr2_diff_list.txt","w")
        print >> diff1_out, "\t".join([str(x['sgr1_diff']) for x in enr_sgr1_sorted])
        print >> diff2_out, "\t".join([str(x['sgr2_diff']) for x in enr_sgr1_sorted])        
        plt.grid(True)
    if plot_type in (1,2,3):
        plotfile = outbase + "_In-subtracted_plot" + plot_type + ".svg"
        plt.savefig(plotfile,format="svg")
        plt.show()

    ### calculate some correlation stats
    
    #(sp_corr,sp_pval) = st.spearmann(x['sgr1_diff'],x['sgr2_diff'])
    #(ks_Dval,ks_pval) = st.ks_2samp(x['sgr1_diff'],x['sgr2_diff'])
    print x['sgr1_diff']
    enrich = [abs(float(x['sgr1_diff']) / x['sgr2_diff']) for x in diff_sgr1_sorted if x['sgr2_diff'] !=0]
    print enrich
    median_enrich = numpy.mean(enrich)
    #print "Spearmann: Corr = %f ; p-val = %f" % (sp_corr,sp_pval)
    #print "KS Test (sp1,sp2): D-statistic = %f ; pval = %f" % (ks_Dval,ks_pval)
    print "Median enrichment: %d ; excluded points: %d" % (median_enrich,len(diff_sgr1_sorted) - len(enrich))

    stat_out = open(outbase + "_sp1-sp2_sgr-diff","w")
    enr_out = open(outbase + "_sp1-sp2_sgr-enr","w")
    
    print >> stat_out, "\t".join([str(x['sgr1_diff']) for x in diff_sgr1_sorted])
    print >> stat_out, "\t".join([str(x['sgr2_diff']) for x in diff_sgr2_sorted])
    
    print >> enr_out, "\t".join([str(x['sgr1_enr']) for x in enr_sgr1_sorted])
    print >> enr_out, "\t".join([str(x['sgr2_enr']) for x in enr_sgr2_sorted])

    

if __name__ == "__main__":
    main(sys.argv[1:])
