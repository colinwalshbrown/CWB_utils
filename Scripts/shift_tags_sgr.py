#!/usr/bin/env python
"""
shift_tags_sgr.py
shift mapped tags (in bowtie format) by a given value, then coallate into an sgr file
"""

import sys

verbose=False

def _main(args):
    
    if len(args) < 4:
        print "usage: shift_tags_sgr.py <shift_size> <output_window> <bowtie_map> <tag_limit (0=none)>--extend"
        sys.exit(0)

    win = int(args[1])
    tags = {}
    chrs = set([])
    tcount = 0
    tag_lim = int(args[3])
    extend = False
    shift = int(args[0])

    if "--extend" in args:
        extend = True
    for line in open(args[2]):
        """
        make a sparse dict of tags
        """
        if tcount % 10000 == 0: 
            print >> sys.stderr, "ntags: %s" % (tcount)
        if tag_lim and (tcount > tag_lim):
            break
        spl = line[:-1].split()
        #print spl
        chr = spl[2]+"_"+spl[1]
        chrs.add(spl[2])
        loc = int(spl[3])
        #strand = spl[1]
        #sh = shift # note tag length not included
        #step = 1
        #if strand == "-":
        #    sh = -shift
        #    step = -1
        if chr in tags:
            #if (extend):
            #    for x in range(loc,loc+sh,step):
            #        #print "%i %i %i" %(loc,sh,step)
            #        addtag(tags,chr,x)
            #else:
            addtag(tags,chr,loc)
        else:
            tags[chr] = {loc:1}
            #if (extend):
            #    for x in range(loc+step,loc+sh,step):
            #        #print "%i %i %i" %(loc,loc+sh,step) 
            #addtag(tags,chr,x)
        #print "%s" % ('\b'* (len(str(tcount)) - 1)),
        tcount += 1

    #print >> sys.stderr, chrs
    #print >> sys.stderr, tags.keys()
    #print >> sys.stderr, tags.items()
    #print >> sys.stderr, sum(map(lambda x: sum(x.values()),tags.values()))
    #print >> sys.stderr, sum(map(sum,tags.values()))
    print >>sys.stderr, tcount
    grand_tot = 0
    for chr in chrs:
        #chr_cnts = []
        plus_chr_name = chr+"_"+"+"
        minus_chr_name = chr+"_"+"-"
        plus_chr={}
        if plus_chr_name in tags:
            plus_chr = tags[plus_chr_name]
        minus_chr={}
        if minus_chr_name in tags:
            minus_chr = tags[minus_chr_name]
#        print "+ chr= ",
#        print plus_chr.keys()
#        print "- chr= ",
#        print minus_chr.keys()
        chr_tot = 0
        chr_max = max(plus_chr.keys() + minus_chr.keys()) + shift
        for x in range(0,chr_max,win):
            tot = 0
            if not extend:
                for p in range(x-shift,(x-shift)+win):
#                    print range(x-shift,(x-shift)+win)
                    if p < 0:
                        continue
                    if p in plus_chr:
#                        print >> sys.stderr, "adding %s:%i to %i => %i" % (plus_chr_name,plus_chr[p],x,tot) 
                        tot += plus_chr[p]
                for m in range(x+(shift-win),x+shift):
#                    print range(x+(shift-win),x+shift)
                    if m in minus_chr:
#                        print >> sys.stderr, "adding %s:%i to %i => %i" % (minus_chr_name,minus_chr[m],x,tot) 
                        tot += minus_chr[m]
            else:
                for p in range((x-shift)+1,x+1):
                    if p < 0:         
                        continue
                    if p in plus_chr:
#                        print >> sys.stderr, "adding %s:%i to %i => %i" % (plus_chr_name,plus_chr[p],x,tot) 
                        tot += plus_chr[p]
                for m in range(x,x+shift):
                    if m in minus_chr:
#                        print >> sys.stderr, "adding %s:%i to %i => %i" % (minus_chr_name,minus_chr[m],x,tot)
                        tot += minus_chr[m]

            print "%s\t%i\t%i" % (chr,x,tot)
            grand_tot += tot
            chr_tot += tot
        print >> sys.stderr, "chr %s total: %i / %i" % (chr,chr_tot,sum(plus_chr.values()) + sum(minus_chr.values()))
    print >> sys.stderr, "grand total: %i" % grand_tot

def addtag(tags,chr,loc):
    if loc in tags[chr]:
        tags[chr][loc] += 1
        if(verbose):
            print "tags @ %s %i: %i" %(chr,loc,tags[chr][loc])
    else:
        tags[chr][loc] = 1
        if(verbose):
            print "tags @ %s %i: %i" %(chr,loc,tags[chr][loc])

if __name__ == "__main__":
    _main(sys.argv[1:])
