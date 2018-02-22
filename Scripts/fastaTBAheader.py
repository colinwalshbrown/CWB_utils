#!/usr/bin/env python

import sys
import re

files = sys.argv[1:]

for fname in files:
    seq = ""
    name = ""
    for line in open(fname):
        headerRE = re.search("^>(\S+)",line)
        if (headerRE):
            if (name):
                newheader = "%s:%s:1:+:%d" % (".".join(fname.split("/")[-1].split(".")[:-1]),name,len(seq))
                print ">" + newheader
                print seq                
            name = headerRE.group(1)
            seq = ""
        elif name:
            seq += line[:-1]
        else:
            continue

