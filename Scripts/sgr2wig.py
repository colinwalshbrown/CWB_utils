#!/usr/bin/env python
##################################
#
# sgr2wig.py convert sgr files to wig (duh!).  Can specify track attributes if needed
#

import sys
from optparse import OptionParser

def _main(args):
    usage = "usage: sgr2wig.py [options] <sgrfile>"
    parser = OptionParser(usage=usage)
    parser.add_option("--name","-n")
    parser.add_option("--viewlimits","-v",help="min:max")
    parser.add_option("--color","-c")
    parser.add_option("--description","-d")
    parser.add_option("--altcolor","-a")
    parser.add_option("--exclude","-e",help="comma-separated list of chr names to exclude")
    parser.add_option("--ucscfix","-u",action="store_true",default=True)
    options,arguments = parser.parse_args()
    if len(arguments) < 1:
        parser.print_help()
        sys.exit(0)

    file = arguments[0]
    name = (".").join(file.split("/")[-1].split(".")[:-1])
    description = ""
    color = "200,100,0"
    altcolor = "100,200,0"
    limits = "0:100"
    exclude = []

    if options.name:
        name = options.name
    if options.description:
        description = options.description
    if options.color:
        color = options.color
    if options.altcolor:
        altcolor = options.altcolor
    if options.exclude:
        exclude = options.exclude.split(",")
    if options.viewlimits:
        limits = options.viewlimits

    trackline = "track type=wiggle_0 name=%s description=%s visibility=full viewLimits=%s color=%s altColor=%s priority=20" % (name,description,limits,color,altcolor)
    print trackline

    prevloc = None
    prevval = None
    prevchr = None
    step = None

    for sgrline in open(file):
        (chr,loc,val) = sgrline[:-1].split()
        if chr in exclude:
            continue
        if prevloc and (chr == prevchr):
            print prevval
        elif not prevloc:
            prevloc = loc
            prevval = val
            continue
        if not (chr == prevchr):
            chrname = chr
            newloc = int(loc)
            if options.ucscfix:
                chrname = "chr" + chr
                newloc += 1
            if not step:
                step = int(loc) - int(prevloc)
            print "fixedStep chrom=%s start=%d step=%s" % (chrname,newloc,step)
        prevloc=loc
        prevchr=chr
        prevval=val

#    for sgrline in open(file):
#        (chr,loc,val) = sgrline[:-1].split()
#        print "%s\t%s\t%s\t%s" % (chr,loc,loc,val)

if __name__ == "__main__":
    _main(sys.argv[1:])

