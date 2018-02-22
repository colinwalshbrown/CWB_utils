#!/usr/bin/env python

import sys

if len(sys.argv) < 2:
    print "usage: TK_xls-to-peaks.py <TK_xls>"
    sys.exit(0)

for line in open(sys.argv[1]):
    l = line[:-1].split()
    for (i,x) in enumerate(l[4][:-1].split(",")):
        print "\t".join((l[0],l[1],x,l[5][:-1].split(",")[i]))
        
