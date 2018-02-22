#!/usr/bin/env python
##################################
#
# sgr2bedGraph.py convert sgr files to begGraph (duh!).  Can specify track attributes if needed
#

import sys
from optparse import OptionParser

def _main(args):
    usage = "usage: sgr2bedGraph.py [options] <sgrfile>"
    parser = OptionParser(usage=usage)
    parser.add_option("--name","-n")
    parser.add_option("--color","-c")
    parser.add_option("--description","-d")
    parser.add_option("--altcolor","-a")
    options,arguments = parser.parse_args()
    if len(arguments) < 1:
        parser.print_help()
        sys.exit(0)

    file = arguments[0]
    name = (".").join(file.split("/")[-1].split(".")[:-1])
    description = ""
    color = "200,100,0"
    altcolor = "100,200,0"

    if options.name:
        name = options.name
    if options.description:
        description = options.description
    if options.color:
        color = options.color
    if options.altcolor:
        altcolor = options.altcolor

    trackline = "track type=bedGraph name=%s description=%s visibility=full color=%s altColor=%s priority=20" % (name,description,color,altcolor)
    print trackline

#     prevloc = None
#     prevval = None
#     prevchr = None

#     for sgrline in open(name):
#         (chr,loc,val) = sgrline[:-1].split()
#         if prevloc and (chr == prevchr):
#             print "%s\t%s\t%s\t%s" % (chr,prevloc,loc,prevval)

    for sgrline in open(file):
        (chr,loc,val) = sgrline[:-1].split()
        print "%s\t%s\t%s\t%s" % (chr,loc,loc,val)

if __name__ == "__main__":
    _main(sys.argv[1:])

