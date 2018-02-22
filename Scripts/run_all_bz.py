#!/usr/bin/env python

import sys
from multiprocessing import Pool
from subprocess import Popen,PIPE

def _main(args):

    if len(args) < 1:
        print "usage: run_all_bz.py <n_processes>"
        sys.exit(1)

    pool = Pool(int(args[0]))
    done = pool.map(build_process,sys.stdin)
    for x in done:
        print x

def build_process(process_string):
    split_out = process_string[:-1].split(" > ")
    out = None
    print split_out
    if len(split_out) > 1:
        out = open(split_out[-1],"w")
    else:
        out = sys.stdout
    cmds = split_out[0].split(" | ")
    print cmds
    p_objs = []
    if len(cmds) > 1:
        for i,c in enumerate(cmds):
            cmd_args = c.split()
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
