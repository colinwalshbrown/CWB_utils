#!/usr/bin/env python

import sys
import sqlite3
import random

def _main(args):
    
    if len(args) < 5:
        print "usage: patser_hit_window_xls.py <patser_hit_db1:patser_hit_db2> <matrix_name1:matrix_name2> <winsize> <xls> <nrandom>"
        sys.exit(1)
    
    dbs = args[0].split(':')
    mtxs = args[1].split(':')
    print mtxs

    (conn,conn2,cur,cur2) = (None,) * 4

    if len(dbs) > 1:
        conn = sqlite3.connect(dbs[0])
        conn2 = sqlite3.connect(dbs[1])
        conn.row_factory = sqlite3.Row
        conn2.row_factory = sqlite3.Row
        cur = conn.cursor()
        cur2 = conn2.cursor()
    else:
        conn = sqlite3.connect(dbs[0])
        conn2 = conn
        conn.row_factory = sqlite3.Row
        cur = conn.cursor()
        cur2 = cur

    winsize = int(args[2])

    xls_name = args[3].split('/')[-1].split(".")[0]
    xls_lines = []
    for line in open(args[3]):
        spl_line = tuple(line[:-1].split("\t"))
        xls_lines.append(spl_line)

    outhits = open(xls_name + "_" + "+".join(mtxs) + "_win" + args[2],"w")
    outrand = open("+".join(mtxs) + "_" + "_randwin" + args[2],"w")

    cur.execute("""SELECT matrix_key from matrix where name = ?""",(mtxs[0],))
    matrix_key = cur.fetchall()[0][0]
    
    cur2.execute("""SELECT matrix_key from matrix where name = ?""",(mtxs[1],))
    matrix2_key = cur2.fetchall()[0][0]

    chrs = {}
    count = 0
    cur.execute("""SELECT DISTINCT chr FROM patser_hit""")
    chr_names = set(cur.fetchall())
    cur2.execute("""SELECT DISTINCT chr FROM patser_hit""")
    chr_names = chr_names.union(cur2.fetchall())
    for ch in chr_names:
        chr = ch[0]
        cur.execute("""SELECT max(start) from patser_hit where chr = ?""",(chr,))
        max = 0
        f = cur.fetchall()[0][0]
        if f:
            max = int(f)
        cur2.execute("""SELECT max(start) from patser_hit where chr = ?""",(chr,))
        f2 = cur2.fetchall()[0][0]
        if (f2) and (int(f2)) > max:
            max = int(f2)
        if (max == 0):
            print chr + " " + str(f) + " " + str(f2)
        chrs[(count,count+max)] = {'name' : chr,
                                   'max' : max }
        count = count + max

#    xls-mt1_counts = [0] * winsize
#    xls-mt2_counts = [0] * winsize
#    xls-mtboth_counts = [0] * winsize
    xls_total = len(xls_lines) 
#    mt1_hit_count = 0
#    mt2_hit_count = 0
    joint_hit_count = 0

    for x in xls_lines:
        peak_loc = int(x[1]) + int(x[4])
        chr = x[0]
        print (chr,peak_loc - (winsize/2),peak_loc + (winsize/2))
        cur.execute("""SELECT * from patser_hit where ((chr = ?) and (start between ? and ?) and (matrix_key = ?))""",(chr,peak_loc - (winsize/2),peak_loc + (winsize/2),matrix_key))
        mtx1_hits = cur.fetchall()
       #map(lambda x: x.update(('mtx'),mtxs[0]), cur.fetchall())
        cur2.execute("""SELECT * from patser_hit where ((chr = ?) and (start between ? and ?) and (matrix_key = ?))""",(chr,peak_loc - (winsize/2),peak_loc + (winsize/2),matrix2_key))
        mtx2_hits = cur2.fetchall()
        for h in (mtx1_hits):
            win_loc = ((peak_loc - h['start']) + winsize / 2) - 1
            hitline = "\t".join(map(str,[h['chr'],h['start'],h['end'],h['strand'],h['score'],h['pval'],mtxs[0],win_loc,'peak',xls_total]))
            print hitline
            print >> outhits, hitline
        for h in (mtx2_hits):
            win_loc = ((peak_loc - h['start']) + winsize / 2) - 1
            hitline = "\t".join(map(str,[h['chr'],h['start'],h['end'],h['strand'],h['score'],h['pval'],mtxs[1],win_loc,'peak',xls_total]))
            print hitline
            print >> outhits, hitline
            
    #rand_counts = [0] * winsize
    #rand_hit_count = 0
    #print xrange(0,10)
    for x in xrange(0,int(args[4])):
        print "random: " + str(x)
        rand = random.randint(0,count)
        chr_rand = None
        sel_chr = None
        range = None
        for c in chrs.keys():
            if rand >= c[0] and rand < c[1]:
                sel_chr = chrs[c]
            else:
                continue
#            temp_count = [0] * winsize
        chr_rand = random.randint(0,sel_chr['max'])
        #print sel_chr
        cur.execute("""SELECT * from patser_hit where ((chr = ?) and (start between ? and ?) and (matrix_key = ?))""",(sel_chr['name'],chr_rand - (winsize/2),chr_rand + (winsize/2),matrix_key))
        mtx1_hits = cur.fetchall()#map(lambda x: x.update(('mtx'),mtxs[0]), cur.fetchall())
        cur2.execute("""SELECT * from patser_hit where ((chr = ?) and (start between ? and ?) and (matrix_key = ?))""",(sel_chr['name'],chr_rand - (winsize/2),chr_rand + (winsize/2),matrix2_key))
        mtx2_hits = cur2.fetchall()#map(lambda x: x.update(('mtx'),mtxs[1]), cur2.fetchall())
        for h in (mtx1_hits):
            win_loc = ((chr_rand - h['start']) + winsize / 2) - 1
            hitline = "\t".join(map(str,[h['chr'],h['start'],h['end'],h['strand'],h['score'],h['pval'],mtxs[0],win_loc,"random",xls_total]))
            #print hitline
            print >> outrand, hitline
        for h in (mtx2_hits):
            win_loc = ((chr_rand - h['start']) + winsize / 2) - 1
            hitline = "\t".join(map(str,[h['chr'],h['start'],h['end'],h['strand'],h['score'],h['pval'],mtxs[1],win_loc,"random",xls_total]))
            #print hitline
            print >> outrand, hitline

#             for h in hits:
#                 rand_hit_count += 1
#                 win_loc = ((chr_rand - h['start']) + winsize / 2) - 1
#                 print win_loc
#                 rand_counts[win_loc] += 1
#                 temp_count[win_loc] += 1
#             print " ".join(map(str,temp_count))
#             print >> outrand, "\t".join(map(str,temp_count))
            
#     countout = open(args[1] + "_win" + args[2] + "counts.txt","w")
#     rateout = open(args[1] + "_win" + args[2] + "rates.txt","w")
#     randrateout = open(args[1] + "_win" + args[2] + "_rep" + args[4] + "rates.txt","w")
#     print "peak windows with at least one hit: %i / %i" % (hit_count,xls_total)
#     print "random windows with at least one hit: %i / %s" % (rand_hit_count,args[4]) 
#     print >> countout, " ".join(map(str,xls_counts))
#     print >> rateout, " ".join([str(x/float(xls_total)) for x in xls_counts])
#     print >> randrateout, " ".join([str(x/float(args[4])) for x in rand_counts])

# def get_joint_hits(mtx1_hits,mtx2_hits,mt1_counts,mt2_counts,joint_counts,loc):
#     temp_count1 = [0] * winsize
#     temp_count2 = [0] * winsize
#     first = False
#     print hits
#     for h1 in mtx1_hits:
#         win_loc1 = ((loc - h1['start']) + winsize / 2) - 1 
#         print win_loc
#         temp_count1[win_loc1] += 1
#         mt1_counts[win_loc1] += 1
#         print " ".join(map(str,temp_count1))
#         #print >> outhits, "\t".join(map(str,temp_count))

#         for h2 in mtx2_hits:
#             win_loc2 = ((loc - h2['start']) + winsize / 2) - 1
#             if not (temp_count2[win_loc2]):
#                 mt2_counts[win_loc2] += 1
#                 temp_count2[win_loc2] += 1
#             print " ".join(map(str,temp_count2))
#             if ((((h1['strand'] == "+") and
#                   (h2['strand'] == "-")) or
#                  ((h1['strand'] == "-") and
#                   (h2['strand'] == "+"))) and
#                 (h1['end'] == h2['start'])):
#                 print "pair:"
#                 print h1
#                 print h2
#                 print "=" * 40
    

            #print >> outhits, "\t".join(map(str,temp_count))


if __name__ == "__main__":
    _main(sys.argv[1:])
