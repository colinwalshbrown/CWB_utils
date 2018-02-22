#!/usr/bin/env python
#
# Wrapper for classes and methods in the fasta_subseq_2.py module; get a single region from a given fasta file
#

import sys
from fasta_subseq_2 import *

def _main(args):

    if len(args) < 4:
        print "usage: fasub.py <fastafile> <chromosome> <start> <end>"
        sys.exit(0)

    start = int(args[2])
    end = int(args[3])

    db = FastaDB()
    db.openFastaFile(args[0])
    sequence = db[args[1]]['sequence'][start:end]
    fa_db_name = args[0].split("/")[-1].split(".")[0]
    print ">%s:%s:%d-%d\n%s" % (fa_db_name,args[1],start,end,sequence)

if __name__ == "__main__":
    _main(sys.argv[1:])
