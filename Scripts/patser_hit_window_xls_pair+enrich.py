#!/usr/bin/env python

import sys
import sqlite3
import random
from subprocess import Popen,PIPE
import fasta_subseq_2
import numpy as np

def _main(args):
    
    if len(args) < 6:
        print "usage: patser_hit_window_xls.py <patser_hit_db1:patser_hit_db2> <matrix_name1:matrix_name2> <sgr1> <fasta> <winsize> <sgr_winsize> <xls>"
        sys.exit(1)
    
    dbs = args[0].split(':')
    mtxs = args[1].split(':')
    sgr = args[2]
    fasta = args[3]
    sgr_winsize = int(args[5])
    print >> sys.stderr, mtxs

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

    winsize = int(args[4])

    xls_name = args[6].split('/')[-1].split(".")[0]
    xls_lines = []
    for line in open(args[6]):
        spl_line = tuple(line[:-1].split("\t"))
        xls_lines.append(spl_line)

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
            print >> sys.stderr, chr + " " + str(f) + " " + str(f2)
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
    
    xls_regions = []

    for x in xls_lines:
        peak_loc = int(x[1]) + int(x[4])
        chr = x[0]
        enrich = x[7]
        #print (chr,peak_loc - (winsize/2),peak_loc + (winsize/2))
        cur.execute("""SELECT * from patser_hit where ((chr = ?) and (start between ? and ?) and (matrix_key = ?))""",(chr,peak_loc - (winsize/2),peak_loc + (winsize/2),matrix_key))
        mtx1_hits = cur.fetchall()
       #map(lambda x: x.update(('mtx'),mtxs[0]), cur.fetchall())
        cur2.execute("""SELECT * from patser_hit where ((chr = ?) and (start between ? and ?) and (matrix_key = ?))""",(chr,peak_loc - (winsize/2),peak_loc + (winsize/2),matrix2_key))
        mtx2_hits = cur2.fetchall()
        xls_ln = {'chr' : chr,
                  'peak' : peak_loc,
                  'enrich' : x[7],
                  'mtx1_hits' : mtx1_hits,
                  'mtx2_hits' : mtx2_hits}
        xls_regions.append(xls_ln)

    xls_regions = get_enrich(xls_regions,sgr,sgr_winsize,fasta)

    for x in xls_regions:
        mtx1_hits = x['mtx1_hits_info']
        mtx2_hits = x['mtx2_hits_info']
        print ("=" * 30) + x['chr'] + ":" + str(x['peak']) + ":enrich=" + x['enrich'] + ("=" * 30)
        for hit in (mtx1_hits):
            h = hit['hit_obj']
            win_loc = ((peak_loc - h['start']) + winsize / 2) - 1
            hitline = "\t".join(map(str,[h['chr'],hit['loc'],h['strand'],h['score'],h['pval'],mtxs[0],win_loc,hit['seq'],hit['enrich_md'],hit['enrich_mn']]))
            print hitline
            #print >> outhits, hitline
        for hit in (mtx2_hits):
            h = hit['hit_obj']
            win_loc = ((peak_loc - h['start']) + winsize / 2) - 1
            hitline = "\t".join(map(str,[h['chr'],hit['loc'],h['strand'],h['score'],h['pval'],mtxs[1],win_loc,hit['seq'],hit['enrich_md'],hit['enrich_mn']]))
            print hitline
            #print >> outhits, hitline
            
    #rand_counts = [0] * winsize
    #rand_hit_count = 0

def get_enrich(xls_regs,sgr,winsize,fasta):

    fasta_db = fasta_subseq_2.FastaDB()
    fasta_db.openFastaFile(fasta)
    for reg in xls_regs:
        
        for x in ('mtx1_hits','mtx2_hits'):
            hit_info=[]
            for h in reg[x]:
                hit = {'hit_obj':h}
                width = abs(h['start'] - h['end'])
                if h['strand'] == "+":
                    seq = fasta_db[h['chr']]['sequence'][h['start']:(h['start']+width)]
                    hit['loc'] = h['start']
                else:
            ### !!!!! CHANGE IF FIX HIT DATABASE!!!!
                    seq = fasta_subseq_2.revcomp(fasta_db[h['chr']]['sequence'][h['end']:(h['end']+width)])
                hit['loc'] = h['end']
                hit['nearest'] = (0,0)
                hit['vals'] = []
                hit['seq'] = seq
                hit_info.append(hit)
            reg[x+'_info'] = hit_info
            
    for y in open(sgr):
        (chr,loc,val) = y.split()
        loc = int(loc)
        val = int(val)
        #print chr
        for x in xls_regs:
            for hit_info in ('mtx1_hits_info','mtx2_hits_info'):
                for d in x[hit_info]:
                #print (loc,target_loc)
                    target_loc = d['loc']
                    if (chr == d['hit_obj']['chr']) and (abs(loc - target_loc) < abs(loc - d['nearest'][0])):
                        d['nearest'] = (loc,val)
                
                    if (chr == d['hit_obj']['chr']) and (abs(loc - target_loc) < (winsize / 2)):
                        d['vals'].append(val)
                        print >> sys.stderr, d
    for x in xls_regs:
        for hit_info in ('mtx1_hits_info','mtx2_hits_info'):
            for h in x[hit_info]:
                h['win_mean'] = np.mean(h['vals'])
                h['win_median'] = np.median(h['vals'])
                h['enrich_md'] = h['nearest'][1] / h['win_median']
                h['enrich_mn'] = h['nearest'][1] / h['win_mean']
                print >> sys.stderr, h

    return xls_regs

if __name__ == "__main__":
    _main(sys.argv[1:])
