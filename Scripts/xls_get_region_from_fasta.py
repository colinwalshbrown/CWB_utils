#!/usr/bin/env python

import sys
import fasta_subseq_2

def _main(args):
    
    if len(args) != 3:
        print ("usage: xls_get_region_from_fasta.py <fasta> <xls> <window>")
        sys.exit(1)

    win = int(args[2])

    fasta = fasta_subseq_2.FastaDB()
    fasta.openFastaFile(args[0])

    seqs = []

    for ln in open(args[1]):
        sp = ln[:-1].split()
        print sp
        pk = int(sp[1]) + int(sp[4])
        seq = fasta[sp[0]]['sequence'][(pk - win):(pk+win)]
        get_in = sp[-1]#raw_input(">%s:%d..%d\n\'k\'=keep; \'r\' = reverse comp; \'<anything else>\' = discard: " % (sp[0],pk-win,pk+win))
        if get_in == 'k':
            pass
        elif get_in == 'r':
            seq = fasta_subseq_2.revcomp(seq)
        else:
            continue

        seqs.append(">%s:%d..%d\n%s" % (sp[0],pk-win,pk+win,seq))

    outfile = raw_input("name of output file: ")
    outfh = open(outfile,"w")
    for s in seqs:
        print >> outfh, s

if __name__ == "__main__":
    _main(sys.argv[1:])
