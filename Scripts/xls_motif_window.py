#!/usr/bin/env python

import sys
import re
import patser_tools
import fasta_subseq_2

def _main(args):

    if len(args) < 4:
        print >> sys.stderr, "usage: xls_motif_window.py <xls> <fasta> <matrix_file> <window>"
        sys.exit(1)

    fasta = fasta_subseq_2.FastaDB()
    fasta.openFastaFile(args[1])
    xls_regions = []
    for x in open(args[0]):
        spl = x[:-1].split()
        region = {'chr':spl[0],'start':int(spl[1]),'end':int(spl[2]),'enrich':spl[7]}
        region['seq'] = fasta[region['chr']]['sequence'][region['start']:region['end']]
        xls_regions.append(region)
        
    for r in xls_regions:
        try:
            annot = patser_tools.makePatserAnnotation(sequence=r['seq'],matrix=args[2])
        except IOError:
            print >>sys.stderr, "Error in seq %s:%d..%d:" % (r['chr'],r['start'],r['end'])
            continue
        if len(annot.getAllFeatures()) < 1:
            continue
        maxhit = annot.getMaxFeature("score")
        winstart = None
        winend = None
        winseq = None
        if maxhit.tags["strand"] == '+':
            winstart = r['start'] + (maxhit.start - int(args[3])/2)
            winend = r['start'] + (maxhit.start + int(args[3])/2)
            win_seq = fasta[r['chr']]['sequence'][winstart:winend]
        else:
            winstart = r['start'] + ((maxhit.end - 3) - int(args[3])/2)
            winend = r['start'] + ((maxhit.end - 3) + int(args[3])/2)
            win_seq = fasta_subseq_2.revcomp(fasta[r['chr']]['sequence'][winstart:winend])
        print ">%s:%d..%d:%s enr=%s mtx=%s" % (r['chr'],winstart,winend,maxhit.tags['strand'],r['enrich'],maxhit.tags['score'])
        print win_seq
        
if __name__=="__main__":
    _main(sys.argv[1:])
