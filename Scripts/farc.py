#!/usr/bin/env python

import sys
import re
from fasta_subseq_2 import revcomp

def _main(args):

    if len(args) < 1:
        print "usage: revcomseq.py <fasta>"
        sys.exit(0)

    sq = ""
    header = ""
    for l in open(args[0]):
        if re.search("^>",l):
            header = l[:-1] + "_rc"
            sq=""
            if (sq):
                tr_seq = re.sub("\n","",sq)
                print "%s\n%s" % (header,revcomp(tr_seq))
        else:
            sq += l
    tr_seq = re.sub("\n","",sq)
    print "%s\n%s" % (header,revcomp(tr_seq))

if __name__ == "__main__":
    _main(sys.argv[1:])
    
