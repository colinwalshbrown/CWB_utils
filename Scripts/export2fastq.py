#!/usr/bin/python
"""
export2fastq.py - convert illumina export sequence files to fastq
"""

import sys

def main(args):
    
    if len(args) < 1:
       print "usage: export2fastq.py <export_file>"
       sys.exit(0)

    for line in open(args[0]):
       line_spl = line[:-1].split()
       fqstr = "%s_%04d:%s:%s:%s:%s#%s/%s" % (line_spl[0],int(line_spl[1]),line_spl[2],line_spl[3],line_spl[4],line_spl[5],line_spl[6],line_spl[7])
       print "@" + fqstr
       print line_spl[8]
       print "+" + fqstr
       print line_spl[9] 

if __name__ == "__main__":
    main(sys.argv[1:])
