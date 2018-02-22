#!/usr/bin/env python

import sys

if len(sys.argv) < 3:
    print "usage: xls_overlap.py <xls1> <xls2> [over_xls_out]"
    sys.exit(0)

xcol = 1;

xls1 = []
xls2 = []
xls_sm = []
xls_lg = []
xls_over = set([])

xls1_over = 0

for line1 in open(sys.argv[1]):
    if line1[0] == "#":
        continue
    spline1 = tuple(line1[:-1].split("\t"))
    xls1.append(spline1)

for line2 in open(sys.argv[2]):
    if line2[0] == "#":
        continue
    spline2 = tuple(line2[:-1].split("\t"))
    xls2.append(spline2)

if len(xls1) > len(xls2):
    xls_sm = xls2
    xls_lg = xls1
else:
    xls_sm = xls1
    xls_lg = xls2

for x1 in xls_sm:
    for x2 in xls_lg:
        if ((((int(x1[1+xcol]) < int(x2[2+xcol])) and (int(x1[1+xcol]) > int(x2[1+xcol]))) or
             ((int(x1[2+xcol]) < int(x2[2+xcol])) and (int(x1[2+xcol]) > int(x2[1+xcol]))) or 
             ((int(x2[1+xcol]) < int(x1[2+xcol])) and (int(x2[1+xcol]) > int(x1[1+xcol]))) or
             ((int(x2[2+xcol]) < int(x1[2+xcol])) and (int(x2[2+xcol]) > int(x2[2+xcol])))) and
            (x1[0+xcol] == x2[0+xcol])):
            xls1_over += 1
            if ((int(x1[2+xcol]) - int(x1[1+xcol])) > (int(x2[2+xcol]) - int(x2[1+xcol]))):
                xls_over.add(x2)
            else:
                xls_over.add(x1)

xls_sm_set = set(tuple(xls_sm))
xls_sm_set.update(xls_lg)
xls_nonover = xls_sm_set.difference(xls_over)

print >> sys.stderr, "xls1: %i of %i total overlap = %f%%" % (xls1_over,len(xls1),((float(xls1_over) / len(xls1)) * 100))
print >> sys.stderr, "xls2: %i of %i total overlap = %f%%" % (xls1_over,len(xls2),(float(xls1_over) / len(xls2)) * 100)

xls_out = sys.stdout

if (len(sys.argv) == 4):
    xls_out = open(sys.argv[3],"w")

for x in xls_over:
    print >> xls_out, "\t".join(x)
xls_out.close()
