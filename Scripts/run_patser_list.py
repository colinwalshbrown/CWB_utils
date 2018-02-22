#!/usr/bin/env python

import sys
import os
import subprocess
import tempfile 

CHUNKSIZE = 10000

def _main(args):
    
    if len(args) < 1:
        print "usage: run_patser_list.py [nreps]"
        sys.exit(1)

    home = ""
    mot_len = 15
    mtx1= "~/Data/her_analysis/dmel_her.txt"
    mtx2= "~/Data/her_analysis/dpse_mtx15_mtx.mtx"

    tempfiles = {}
    for x in range(nreps/Chunksize):
        (jobid,scriptfile,outfile) = makeJob(mot_len,CHUNKSIZE,home,mtx1,mtx2)
        tempfiles[jobid] = {'script':scriptfile,
                            'out':outfile,
                            'done':False}

def makeJob(mot_len,reps,home,mtx1,mtx2):
    scriptfile = tempfile.NamedTemporaryFile(mode="w+b",delete=False) 
    outfile = tempfile.NamedTemporaryFile(mode="w+b",delete=False)
    script_text = """!#/bin/sh
%s/scripts/random_seqs.py %d %d | %s/scripts/patser_list.py %s %s > %s""" % (home,mot_len,reps,home,mtx1,mtx2,outfile.name)
    print >> scriptfile, script_text
    scriptfile.close()
    jobout = Popen(xgrid
    return (scriptfile,outfile)

if __name__ == "__main__":
    _main(sys.argv[1:])
