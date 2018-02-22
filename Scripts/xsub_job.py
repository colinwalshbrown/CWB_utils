#!/usr/bin/python

import sys
import os
import re
import time
import xgrid_tools_2
from subprocess import Popen,PIPE

CHK_INTERVAL = 1

def _main(args):

    dir = None
    verbose = False
    cmds = []
    inter = 0

    if len(args) < 1:
        print "usage: xsub_job.py (--dir=RUNDIR --outdir=OUTDIR --verbose) [CMD]"
        sys.exit(1)

    for a in args:
        dirsearch = re.search("(--dir|-d)\s+(\S+)",a)
        verbsearch = re.search("(--verbose|-v)",a)
        outsearch = re.search("(--outdir|-o)\s+(\S+)",a)

        if dirsearch:
            dir = dirsearch.groups()[1]
        else:
            dir = os.getcwd()
        
        if outsearch:
            out = open(outsearch.groups()[1],"w")
        else:
            out = sys.stdout

        if verbsearch:
            verbose = True
            
        if not (verbsearch or dirsearch):
            cmds.append(a)
    
    cmd = " ".join(cmds)

    job = xgrid_tools_2.xgrid_job(run_cmd=cmd,dir=dir)
    try:
        job.start()
        if (verbose):
            print >> sys.stderr,"Started job: jobID=" + job.jobID
            
        while job.checkstatus() == "RUNNING":
            time.sleep(CHK_INTERVAL)
            inter += 1
            if (verbose):
                print "Job %s Running %d:%d:%d" % (job.jobID,int(inter/3600),int((inter % 3600)/60),int((inter % 60))) 
            
        if job.checkstatus() == "FAIL":
                print >> sys.stderr,"Job %s Failed" % (str(job.jobID),)
        elif job.checkstatus() == "DONE":
            print >> sys.stderr,"Job %s Done" % (str(job.jobID),)
            print >> out, job.results()
        else:
            print >> sys.stderr, "Job %s Fucked Up Incomprehensibly"
    finally:    
        job.cleanup()


if __name__ == "__main__":
    _main(sys.argv[1:])
