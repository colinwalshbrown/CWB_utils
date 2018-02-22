#!/usr/bin/env python
"""
take file of paired peaks (format: sp1_name[tab]sp1_chr[tab]sp1_peak[tab]sp2_name[tab]sp2chr[tab]sp2_peak); extract alignments window around ea. peak,
search for motif occurances
"""

import sys
import re
import patser_tools
import alignment
import sqlite3

def main(args):
    
    if len(args) < 7:
        print "usage: patser_annotate_paired_peaks.py <paired_peak_file> <alignment> <sp1_patser_db> <sp2_patser_db> <sp1_matrix_name> <sp2_matrix_name> <window> <basename>"
        sys.exit(1)
        
    aln = alignment.WholeGenomeAlign(args[1],format="FSA",reindex=False)
    conn_sp1 = sqlite3.connect(args[2])
    conn_sp1.row_factory = sqlite3.Row
    conn_sp2 = sqlite3.connect(args[3])
    conn_sp2.row_factory = sqlite3.Row
    cur_sp1 = conn_sp1.cursor()
    cur_sp2 = conn_sp2.cursor()
    winsize = int(args[6])
    outbase = args[7]

    mtx1 = args[4]
    mtx2 = args[5]
    peaks = args[0]

    cur_sp1.execute("""SELECT matrix_key from matrix where name = ?""",(mtx1,))
    matrix_key_sp1 = cur_sp1.fetchall()[0][0]

    cur_sp2.execute("""SELECT matrix_key from matrix where name = ?""",(mtx2,))
    matrix_key_sp2 = cur_sp2.fetchall()[0][0]

    print matrix_key_sp1
    print matrix_key_sp2

    sp1_alns_out = open(outbase + "_sp1_alns.txt","w")
    sp2_alns_out = open(outbase + "_sp2_alns.txt","w")
    hits_out = open(outbase + "_pair-hits.txt","w")

    pks = []

    counts = {'sp1_only_win':0,
              'sp1_sp2_win':0,
              'sp2_only_win':0,
              'neither_win':0,
              'sp1_only_site':0,
              'sp1_sp2_site':0,
              'sp2_only_win':0,
              'neither_site':0}

    for ln in open(peaks):
        if re.search("#",ln):
            continue
        ln_spl = ln[:-1].split("\t")
        #pk_hts = ln_spl[5][:-1].split(",")
        #print pk_hts
        for pk in ln_spl[2].split(","):
            if not pk:
                continue
            for pk2 in ln_spl[5].split(","):
                new_ln = [ln_spl[0],ln_spl[1],pk,ln_spl[3],ln_spl[4],pk2]
                pks.append(new_ln)

    for pk in pks:
        sp1_chr = pk[1]
        sp1_loc = int(pk[2])
        sp1_winstart = sp1_loc - (winsize/2)
        sp1_winend = sp1_loc + (winsize/2)
        cur_sp1.execute("""SELECT * from patser_hit where ((chr = ?) and (start between ? and ?) and (matrix_key = ?))""",(sp1_chr,sp1_winstart,sp1_winend,matrix_key_sp1))
        sp1_sp1_hits = cur_sp1.fetchall()
        sp1_sp2_aln = aln.getSliceSpeciesCoord(pk[0],sp1_chr,sp1_winstart,sp1_winend)[0]
        sp1_sp2_chr = sp1_sp2_aln[pk[3]].parent
        sp1_sp2_winstart = sp1_sp2_aln[pk[3]].parentStart
        sp1_sp2_winend = sp1_sp2_aln[pk[3]].parentEnd
        sp1_sp2_strand = sp1_sp2_aln[pk[3]].parentStrand
        cur_sp2.execute("""SELECT * from patser_hit where ((chr = ?) and (start between ? and ?) and (matrix_key = ?))""",(sp1_sp2_chr,sp1_sp2_winstart,sp1_sp2_winend,matrix_key_sp2))
        sp1_sp2_hits = cur_sp2.fetchall()

        
        sp2_chr = pk[4]
        sp2_loc = int(pk[5])
        sp2_winstart = sp2_loc - (winsize/2)
        sp2_winend = sp2_loc + (winsize/2)
        cur_sp2.execute("""SELECT * from patser_hit where ((chr = ?) and (start between ? and ?) and (matrix_key = ?))""",(sp2_chr,sp2_winstart,sp2_winend,matrix_key_sp2))
        sp2_sp2_hits = cur_sp2.fetchall()
        sp2_sp1_aln = aln.getSliceSpeciesCoord(pk[3],sp2_chr,sp2_winstart,sp2_winend)[0]
        sp2_sp1_chr = sp2_sp1_aln[pk[0]].parent
        sp2_sp1_winstart = sp2_sp1_aln[pk[0]].parentStart
        sp2_sp1_winend = sp2_sp1_aln[pk[0]].parentEnd
        sp2_sp1_strand = sp2_sp1_aln[pk[0]].parentStrand
        print("""SELECT * from patser_hit where ((chr = %s) and (start between %s and %s) and (matrix_key = %s))""" %(sp2_sp1_chr,sp2_sp1_winstart,sp2_sp1_winend,matrix_key_sp1))

        cur_sp1.execute("""SELECT * from patser_hit where ((chr = ?) and (start between ? and ?) and (matrix_key = ?))""",(sp2_sp1_chr,sp2_sp1_winstart,sp2_sp1_winend,matrix_key_sp1))
        sp2_sp1_hits = cur_sp1.fetchall()
        
        
        sp1_sp1_str = map(lambda x: "%s-%s:%d-%d" % (pk[0],x['chr'],x['start'],x['end']),sp1_sp1_hits)
        if (sp1_sp2_winstart < sp1_sp2_winend):
            sp1_sp2_str = map(lambda x: "%s-%s:%d-%d:%s" % (pk[3],sp1_sp2_aln[pk[0]].parent,sp1_sp2_aln.transformSpCoord(pk[3],pk[0],x['start']),sp1_sp2_aln.transformSpCoord(pk[3],pk[0],x['end']),sp1_sp2_strand),sp1_sp2_hits)
        else:
            sp2_sp1_str = map(lambda x: "%s-%s:%d-%d" % (pk[0]+"(no orth seq)",x['chr'],x['start'],x['end']),sp1_sp2_hits)
        sp2_sp1_str = map(lambda x: "%s-%s:%d-%d" % (pk[0],x['chr'],x['start'],x['end']),sp2_sp1_hits)
        print sp1_sp2_aln
        print sp2_sp1_aln
        if (sp2_sp1_winstart < sp2_sp1_winend):
            sp2_sp2_str = map(lambda x: "%s-%s:%d-%d:%s" % (pk[3],sp2_sp1_aln[pk[0]].parent,sp2_sp1_aln.transformSpCoord(pk[3],pk[0],x['start']),sp2_sp1_aln.transformSpCoord(pk[3],pk[0],x['end']),sp1_sp2_strand),sp2_sp2_hits)
#            print "FAIL: "
        else:
            sp2_sp2_str = map(lambda x: "%s-%s:%d-%d" % (pk[3] + "(no orth seq)",x['chr'],x['start'],x['end']),sp2_sp2_hits)
        print sp1_sp1_hits
        print sp1_sp2_hits
        print sp2_sp1_hits
        print sp2_sp2_hits

        print >> sp1_alns_out, sp1_sp2_aln
        print >> sp2_alns_out, sp2_sp1_aln

#        sp1_sp2_aln.PrintClustal(100,outHandle=sp1_alns_out)
#        sp2_sp1_aln.PrintClustal(100,outHandle=sp2_alns_out)
                          
        print >> hits_out, "%s:%s:%s-%s %s %s %s:%s:%s-%s %s %s" % (pk[0],sp1_chr,sp1_winstart,sp1_winend,sp1_sp1_str,sp1_sp2_str,pk[3],sp2_chr,sp2_winstart,sp2_winend,sp2_sp2_str,sp2_sp1_str)
        #print sp1_sp2_aln

if __name__ == "__main__":
    main(sys.argv[1:])
