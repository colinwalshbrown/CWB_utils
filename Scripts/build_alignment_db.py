#!/usr/bin/env python

import sys
import alignment

def _main(args):

    if len(args) != 1:
        print "USAGE: build_alignment_db.py <fsa_aln_dir_location>"
        sys.exit(0)
    print args[0]
    new_aln = alignment.WholeGenomeAlign(alnFile=args[0] + "/map",format="fsa",reindex=True)
    

if __name__ == "__main__":
    _main(sys.argv[1:])
