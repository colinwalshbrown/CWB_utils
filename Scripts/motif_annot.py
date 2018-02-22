#!/usr/bin/env python

import sys
import re
import patser_tools

def _main(args):
    
    if len(args) < 1:
        print >> sys.stderr, "usage: <mtx_1> [<mtx2>...] < fasta"
        sys.exit(1)
        
    seqname = ""
    seq = ""
    for x in sys.stdin:
        head_re = re.search(">(\S+)",x)
        if head_re and not seq:
            seqname = head_re.group(1)
            continue
        elif head_re:
            for x in args:
                annot = patser_tools.makePatserAnnotation(sequence=seq,matrix=x)
                for x in annot.getAllFeatures():
                    if x.tags['strand'] == '-':
                        print "%s\t%s\t%d\t%s\t%f" % (seqname,x.tags['motif_name'],x.end,x.tags['strand'],x.tags['score'])
                    elif x.tags['strand'] == '+':
                        print "%s\t%s\t%d\t%s\t%f" % (seqname,x.tags['motif_name'],x.start,x.tags['strand'],x.tags['score'])
                
            seqname = head_re.group(1)
            seq = ""
        else:
            seq += x.rstrip()
    

if __name__ == "__main__":
    _main(sys.argv[1:])
