#!/usr/bin/env python

import sys
import os
import re
from multiprocessing import Pool
from subprocess import Popen,PIPE

XGRID_HOST="kermit.QB3.Berkeley.EDU"
XGRID_RUN= ["xgrid",
            "-h",
            XGRID_HOST,
            "-auth",
            "Kerberos",
            "-job",
            "run",
            "-in",
            os.getcwd()]

def _main(args):

    if len(args) < 1:
        print "usage: run_all_bz.py [input from stdin]"
        sys.exit(1)
    bzw_jobs = [x for x in sys.stdin if re.search("blastzWrapper",x)]
    sng_cov_jobs = [x for x in sys.stdin if re.search("single_cov",x)]
    print bzw_jobs
    print sng_cov_jobs
    pool = Pool(int(args[0]))
    bzw_done = pool.map(build_process,bzw_jobs)
    for x in bzw_done:
        print x + "=> DONE"
    sngc_done = pool.map(build_process,sng_cov_jobs)
    for x in sngc_done:
        print x + "=> DONE"

def build_process(process_string):
    split_out = process_string[:-1].split(" > ")
    out = None
    print split_out
    if len(split_out) > 1:
        out = open(split_out[-1],"w")
    else:
        out = sys.stdout
    cmds = split_out[0].split(" | ")
    #print cmds
    p_objs = []
    if len(cmds) > 1:
        for i,c in enumerate(cmds):
            cmd_args = XGRID_RUN + c.split()
            print cmd_args
            p = None
            if i == 0:
                p = Popen(cmd_args,stdout=PIPE)
            elif i == len(cmds)-1:
                p = Popen(cmd_args,stdin=p_objs[-1].stdout,stdout=out) 
            else:
                p = Popen(cmd_args,stdin=p_objs[-1].stdout,stdout=PIPE)
            print(p,cmd_args)    
            p_objs.append(p)
    else:
        p_objs.append(Popen(cmds[0].split(),stdout=out))
    
    p_objs[-1].communicate()
    done_str = "%s Done" % process_string
    return done_str            

if __name__ == "__main__":
    _main(sys.argv[1:])
