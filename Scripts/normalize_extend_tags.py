#!/usr/bin/env python
#
# normalize_tags.py <norm. factor, default= 1000000> <shifted_sgr_files...>
#
# count total tags in (shifted) sgr file, normalize tags at each position by factor (norm. factor) / total

import sys

def _main(args):

    if len(args) < 2:
        print "usage:\n\nnormalize_shifted_tags.py <norm. factor> <frag_size> <step_size> <shifted_sgr_1> [...]"
        sys.exit(1)

    nfact = args[0]
    frag_size = int(args[1])/2
    step_size = int(args[2])
    for file in args[3:]:
        ntags = 0

        tags = {}
        
        # count tags in file
        for x in open(file):
            line = x[:-1].split()
            ntags += int(line[2])
            if line[0] in tags.keys():
                if tags[line[0]]:
                    tags[line[0]].append(line[2])
                else:
                    tags[line[0]] = [line[2]]
            else:
                print line[0]
                tags[line[0]] = [line[2]]

        print '%s: %d tags' % (file,ntags)

        # normalize tag count at each window
        norm = float(nfact) / ntags
        output = open(file + "_norm" + nfact,'w')
        ext = frag_size / step_size

        normtags = {}

        for (chr,vals) in tags.items():
            for (loc,val) in enumerate(vals):
                #print (loc,val)
                for x in range((loc - ext),(loc + ext)):
                    if x < 0:
                        continue
                    else:
                        newloc = x*step_size
                        if chr in normtags.keys():
                            if newloc in normtags[chr].keys():
                                normtags[chr][newloc] += int(val) * norm
                            else:
                                normtags[chr][newloc] = (int(val)*norm)
                        else:
                            print "extending %s" % (chr,)
                            normtags[chr] = {newloc : (int(val)*norm)}
                        #print (chr,newloc,normtags[chr][newloc])
        print "Tag calc. done; writing"

        for chr in sorted(normtags.keys()):
            print chr
            print len(normtags[chr])
            for (loc,val) in sorted(normtags[chr].items()):
                print >> output, "%s\t%s\t%f" % (chr,loc,normtags[chr][loc])

        output.close()

if __name__ == "__main__":
    _main(sys.argv[1:])
