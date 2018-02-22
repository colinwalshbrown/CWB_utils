#!/usr/bin/env python
"""
Build a fastq file that is a random mix of reads from two input fastq files, with the proportion specified
"""

import sys
from random import random

def main(args):
    
    if len(args) != 4:
        print "usage:\n\nmake_mixed_fastq.py <fastq_file_1> <fastq_file_2> <proportion_fastq_1> <n_output_reads>"
        sys.exit(0)

    fq1 = open(args[0])
    fq2 = open(args[1])

    q = float(args[2]) # proportion of fq1
    n = int(args[3]) # total output sample size

    for n in range(0,n):
        r = random()
        if r <= q:
            print read_n_lines(4,fq1)[:-1]
        else:
            print read_n_lines(4,fq2)[:-1]

def read_n_lines(n,file):
    lines = ""
    for x in range(0,n):
        lines += file.readline()
    return lines

if __name__ == "__main__":
    main(sys.argv[1:])
