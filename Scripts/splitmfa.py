#!/usr/bin/python

import sys
import re

def _main(args):
    cur_file = None
    for line in open(args[0]):
        head_match = re.search(">\s?(\S+)",line)
        if head_match:
            name = head_match.group(1)
            print name + ".fa"
            if cur_file:
                cur_file.close()
            cur_file = open(name + ".fa","w")
        print >> cur_file, line,


if (__name__ == "__main__"):
    _main(sys.argv[1:])
