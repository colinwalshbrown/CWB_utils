#!/usr/bin/env python
"""
symlink-dir.py take a directory as input and one as output; symlink all files in dir1 in dir2 
"""

import os
import os.path
import sys
from itertools import imap

def main(args):
    
    if len(args) < 2:
        print "usage: symlink-dir.py <targets_dir> <links_dir>"
        sys.exit(0)

    targets_dir_path = os.path.abspath(args[0])
    links_dir_path = os.path.abspath(args[1]) 
    
    targets = os.listdir(targets_dir_path)
    target_paths = [os.path.join(targets_dir_path,x) for x in targets]
    link_paths = [os.path.join(links_dir_path,x) for x in targets]
    print target_paths
    print link_paths
    results = map(os.symlink,target_paths,link_paths)
    print results

if __name__ == "__main__":
    main(sys.argv[1:])
