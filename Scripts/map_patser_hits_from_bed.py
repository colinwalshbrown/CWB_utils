#!/usr/bin/env python
"""
map_patser_hits_from_bed.py - 
"""

import sys
import patser_tools
import fasta_subseq_2
from pprint import pprint as pp

def _main(args):
    
    if len(args) != 3:
        print "usage: <bed_file> <seq_file> <matrix>"
        sys.exit(0)

    fasta = fasta_subseq_2.FastaDB()
    fasta.openFastaFile(args[1])

    bed_annots = []
    bed_in = open(args[0])
    
    for line in bed_in:
        
        spl = line[:-1].split()
        fseq = fasta[spl[0]]["sequence"][int(spl[1]):int(spl[2])]
        if spl[5] == "-":
            fseq = fasta_subseq_2.revcomp(fseq)
        #print spl
        try:
            patannot = patser_tools.makePatserAnnotation(sequence=fseq,matrix=args[2])
        except:
            continue
        #print "-" * 30
        #print spl
        #print pp(patannot.getAllFeatures())
        bed_annots.append({"seq" : spl[0] + "_" + spl[1] + "_" + spl[2],
                           "annotation":patannot})
        
    for ann in bed_annots:
        for feat in ann["annotation"].getAllFeatures():
            print "%s\t%i\t%i\t%f\t%f\t%s" % (ann["seq"],feat.st,feat.en,feat.tagset["score"],feat.tagset["pval"],feat.tagset["strand"])
            

if __name__ == "__main__":
    _main(sys.argv[1:])
