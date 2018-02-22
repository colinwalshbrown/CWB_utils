#!/usr/bin/env python

import sys 
import re

def main(args):

    if len(args) < 2:
        print "filter_Tommy_xls_by_MACS_xls.py <MACS_xls> <Tommy_xls>"
        sys.exit(1)

    macs = []
    tommy = []

    for x in open(args[0]):
        if re.search("#",x):
            continue
        spl = x.split()
        macs.append(spl)
    
    for x in open(args[1]):
        if re.search("#",x):
            print "?"
            continue
        spl = x.split()
        tommy.append(spl)

    for x in tommy:
        t_chr = x[1]
        t_start = int(x[2])
        t_end = int(x[3])
        for y in macs:
            m_chr = y[0] 
            m_start = int(y[1])
            m_end = int(y[2])
            
            if ((t_chr == m_chr) and 
                ((((t_start > m_start) and
                   (t_start < m_end)) or
                  ((t_end > m_start) and
                   (t_end < m_start))) or
                 (((m_start > t_start) and
                   (m_start < t_end)) or
                  ((m_end > t_start) and
                   (m_start < t_end))))):
                print "\t".join(x)

if __name__ == "__main__":
    main(sys.argv[1:])
