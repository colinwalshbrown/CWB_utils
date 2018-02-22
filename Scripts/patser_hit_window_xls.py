#!/usr/bin/env python

import sys
import sqlite3
import random

def _main(args):
    
    if len(args) < 6:
        print "usage: patser_hit_window_xls.py <patser_hit_db> <matrix_name> <winsize> <xls> <nrandom> <xls-type; 0=(name\tchr\tpeak);1=MACS>"
        sys.exit(1)

    conn = sqlite3.connect(args[0])
    conn.row_factory = sqlite3.Row
    cur = conn.cursor()

    xls_type = int(args[5])
    winsize = int(args[2])

    xls_lines = []
    for line in open(args[3]):
        spl_line = tuple(line[:-1].split("\t"))
        xls_lines.append(spl_line)

    outhits = open(args[3] + "_" + args[1] + "_win" + args[2],"w")
    outrand = open(args[1] + "_" + "_randwin" + args[2],"w")

    cur.execute("""SELECT matrix_key from matrix where name = ?""",(args[1],))
    matrix_key = cur.fetchall()[0][0]

    chrs = {}
    count = 0
    cur.execute("""SELECT DISTINCT chr FROM patser_hit""")
    chr_names = cur.fetchall()
    for ch in chr_names:
        chr = ch[0]
        cur.execute("""SELECT max(start) from patser_hit where chr = ?""",(chr,))
        max = int(cur.fetchall()[0][0])
        chrs[(count,count+max)] = {'name' : chr,
                                   'max' : max }
        count = count + max

    xls_counts = [0] * winsize
    xls_total = len(xls_lines)
    hit_count = 0
    for x in xls_lines:
        peak_loc = None
        chr = None
        if (xls_type == 1):
            peak_loc = int(x[1]) + int(x[4])
            chr = x[0]
        else:
            peak_loc = int(x[2])
            chr = x[1]
        temp_count = [0] * winsize
        print (chr,peak_loc - (winsize/2),peak_loc + (winsize/2))
        cur.execute("""SELECT * from patser_hit where ((chr = ?) and (start between ? and ?) and (matrix_key = ?))""",(chr,peak_loc - (winsize/2),peak_loc + (winsize/2),matrix_key))
        hits = cur.fetchall()
        print hits
        for h in hits:
            hit_count += 1
            win_loc = ((peak_loc - h['start']) + winsize / 2) - 1 
            print win_loc
            xls_counts[win_loc] += 1
            temp_count[win_loc] += 1
        print " ".join(map(str,temp_count))
        print >> outhits, "\t".join(map(str,temp_count))
    
    rand_counts = [0] * winsize
    rand_hit_count = 0
    print xrange(0,10)
    for x in xrange(0,int(args[4])):
        print "random: " + str(x)
        rand = random.randint(0,count)
        sel_chr = None
        range = None
        for c in chrs.keys():
            if rand >= c[0] and rand < c[1]:
                sel_chr = chrs[c]
            else:
                continue
            temp_count = [0] * winsize
            chr_rand = random.randint(0,sel_chr['max'])
            cur.execute("""SELECT * from patser_hit where ((chr = ?) and (start between ? and ?) and (matrix_key = ?))""",(chr,chr_rand - (winsize/2),chr_rand + (winsize/2),matrix_key))
            hits = cur.fetchall()
            for h in hits:
                rand_hit_count += 1
                win_loc = ((chr_rand - h['start']) + winsize / 2) - 1
                print win_loc
                rand_counts[win_loc] += 1
                temp_count[win_loc] += 1
            print " ".join(map(str,temp_count))
            print >> outrand, "\t".join(map(str,temp_count))
            
    countout = open(args[1] + "_win" + args[2] + "counts.txt","w")
    rateout = open(args[1] + "_win" + args[2] + "rates.txt","w")
    randrateout = open(args[1] + "_win" + args[2] + "_rep" + args[4] + "rates.txt","w")
    print "peak windows with at least one hit: %i / %i" % (hit_count,xls_total)
    print "random windows with at least one hit: %i / %s" % (rand_hit_count,args[4]) 
    print >> countout, " ".join(map(str,xls_counts))
    print >> rateout, " ".join([str(x/float(xls_total)) for x in xls_counts])
    print >> randrateout, " ".join([str(x/float(args[4])) for x in rand_counts])

if __name__ == "__main__":
    _main(sys.argv[1:])
