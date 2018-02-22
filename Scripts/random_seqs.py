#!/usr/bin/env python

import sys
import random

def _main(args):
    
    if len(args) < 2:
        print "usage: random_seqs.py <seq_length> <n_reps>"
        sys.exit(1)

    bases = ["A","C","G","T"]
    for x in range(0,int(args[1])):
        samp = "".join([random.choice(bases) for y in range(int(args[0]))])
        print ">rand%s\n%s" % (x,samp)

    
if __name__ == "__main__":
    _main(sys.argv[1:])
