#!/usr/bin/env python

#
# normalize_tags.py <norm. factor, default= 1000000> <shifted_sgr_files...>
#
# count total tags in (shifted) sgr file, normalize tags at each position by factor (norm. factor) / total

import sys

def _main(args):

    if len(args) < 2:
        print "usage:\n\nnormalize_tags.py <norm_factor> <shifted_sgr_1> [...]"
        sys.exit(1)

    nfact = args[0]
    ntags = 0

    for file in args[1:]:

        # count tags in file
        for x in open(file):
            line = x[:-1].split()
            ntags += int(line[2])

        # normalize tag count at each window
        norm = float(nfact) / ntags
        output = open(file + "_norm" + nfact,'w')
        normtags = 0.0

        for x in open(file):
            line = x[:-1].split()
            nwin = int(line[2]) * norm
            print >> output, "%s\t%s\t%f" % (line[0],line[1],nwin)
            normtags += nwin

        print '%s: %f tag equivalents after normalization' % (file,normtags)

        output.close()

if __name__ == "__main__":
    _main(sys.argv[1:])
