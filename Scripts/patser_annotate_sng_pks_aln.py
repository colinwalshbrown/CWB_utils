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
    
    if len(args) < 10:
        print "usage: patser_annotate_paired_peaks.py <xls> <species1> <species2> <alignment> <sp1_patser_db> <sp1_matrix_name> <sp2_patser_db> <sp2_matrix_name> <window> <basename>"
        sys.exit(1)
        
    print >> sys.stderr, args
    aln = alignment.WholeGenomeAlign(args[3],format="FSA",reindex=False)
    conn_sp1 = sqlite3.connect(args[4])
    conn_sp1.row_factory = sqlite3.Row
    conn_sp2 = sqlite3.connect(args[6])
    conn_sp2.row_factory = sqlite3.Row
    cur_sp1 = conn_sp1.cursor()
    cur_sp2 = conn_sp2.cursor()
    winsize = int(args[8])
    basename = args[9]
    sp1 = args[1]
    sp2 = args[2]

    mtx1 = args[5]
    mtx2 = args[7]
    peaks = args[0]

    cur_sp1.execute("""SELECT matrix_key from matrix where name = ?""",(mtx1,))
    matrix_key_sp1 = cur_sp1.fetchall()[0][0]

    cur_sp2.execute("""SELECT matrix_key from matrix where name = ?""",(mtx2,))
    matrix_key_sp2 = cur_sp2.fetchall()[0][0]

    hits_out = open(basename + str(winsize) + "_hits.txt","w")
    alns_out = open(basename + str(winsize) + "_alns.txt","w")
    counts_out = open(basename + str(winsize) + "_counts.txt","w")

    print matrix_key_sp1
    print matrix_key_sp2

    pks = []

    counts = {'sp1_only_win':0,
              'sp1_sp2_win':0,
              'sp2_only_win':0,
              'neither_win':0,
              'sp1_only_site':0,
              'sp1_sp2_site':0,
              'sp2_only_site':0,
              'neither_site':0}

    for ln in open(peaks):
        if re.search("#",ln):
            continue
        ln_spl = ln[:-1].split(" ")
        #pk_hts = ln_spl[5][:-1].split(",")
        #print pk_hts
        for (idx,pk) in enumerate(ln_spl[4].split(",")):
            if not pk:
                continue
            pks.append({'chr':ln_spl[1],'reg_id':ln_spl[0],'loc':int(pk)})
            print({'chr':ln_spl[1],'reg_id':ln_spl[0],'loc':int(pk)})

    for pk in pks:
        sp1_chr = pk['chr']
        sp1_loc = pk['loc']
        sp1_winstart = sp1_loc - (winsize/2)
        sp1_winend = sp1_loc + (winsize/2)
        cur_sp1.execute("""SELECT * from patser_hit where ((chr = ?) and (start between ? and ?) and (matrix_key = ?))""",(sp1_chr,sp1_winstart,sp1_winend,matrix_key_sp1))
        sp1_sp1_hits = cur_sp1.fetchall()
        sp1_sp2_aln_in = aln.getSliceSpeciesCoord(sp1,sp1_chr,sp1_winstart,sp1_winend)
        if len(sp1_sp2_aln_in) != 1:
            continue
        elif len(sp1_sp2_aln_in) > 1:
            print "split aln: " + str(pk)
        sp1_sp2_aln = sp1_sp2_aln_in[0]
        sp1_sp2_chr = sp1_sp2_aln[sp2].parent
        sp1_sp2_winstart = sp1_sp2_aln[sp2].parentStart
        sp1_sp2_winend = sp1_sp2_aln[sp2].parentEnd
        cur_sp2.execute("""SELECT * from patser_hit where ((chr = ?) and (start between ? and ?) and (matrix_key = ?))""",(sp1_sp2_chr,sp1_sp2_winstart,sp1_sp2_winend,matrix_key_sp2))
        sp1_sp2_hits = cur_sp2.fetchall()

        
#        sp2_chr = pk[4]
#        sp2_loc = int(pk[5])
#        sp2_winstart = sp2_loc - (winsize/2)
#        sp2_winend = sp2_loc + (winsize/2)
#        cur_sp2.execute("""SELECT * from patser_hit where ((chr = ?) and (start between ? and ?) and (matrix_key = ?))""",(sp2_chr,sp2_winstart,sp2_winend,matrix_key_sp2))
#        sp2_sp2_hits = cur_sp2.fetchall()
#        sp2_sp1_aln = aln.getSliceSpeciesCoord(pk[0],sp1_chr,sp1_winstart,sp1_winend)[0]
#        sp2_sp1_chr = sp2_sp1_aln[pk[0]].parent
#        sp2_sp1_winstart = sp2_sp1_aln[pk[0]].parentStart
#        sp2_sp1_winend = sp2_sp1_aln[pk[0]].parentEnd
#        print("""SELECT * from patser_hit where ((chr = %s) and (start between %s and %s) and (matrix_key = %s))""" %(sp2_sp1_chr,sp2_sp1_winstart,sp2_sp1_winend,matrix_key_sp1))

#        cur_sp1.execute("""SELECT * from patser_hit where ((chr = ?) and (start between ? and ?) and (matrix_key = ?))""",(sp2_sp1_chr,sp2_sp1_winstart,sp2_sp1_winend,matrix_key_sp1))
#        sp2_sp1_hits = cur_sp1.fetchall()
        
        sp1_sp1_str = map(lambda x: "%s:%d-%d" % (x['chr'],x['start'],x['end']),sp1_sp1_hits)
        sp1_sp2_str = map(lambda x: "%s:%d-%d" % (x['chr'],x['start'],x['end']),sp1_sp2_hits)
#        sp2_sp2_str = map(lambda x: "%s:%d-%d" % (x['chr'],x['start'],x['end']),sp2_sp2_hits)
#        sp2_sp1_str = map(lambda x: "%s:%d-%d" % (x['chr'],x['start'],x['end']),sp2_sp1_hits)

        if (len(sp1_sp1_str) > 0) and (len(sp1_sp2_str) > 0):
            print 'both'
            counts['sp1_sp2_win'] += 1
            (shared,sp1_only,sp2_only) = comp_site_aln(sp1_sp1_hits,sp1,sp1_sp2_hits,sp2,sp1_sp2_aln)
            counts['sp1_sp2_site'] += shared
            counts['sp1_only_site'] += sp1_only
            counts['sp2_only_site'] += sp2_only
        elif (len(sp1_sp1_str) > 0) and (len(sp1_sp2_str) < 1):
            print sp1
            counts['sp1_only_win'] += 1
        elif (len(sp1_sp1_str) < 1) and (len(sp1_sp2_str) > 0):
            print sp2
            counts['sp2_only_win'] += 1
        elif (len(sp1_sp1_str) == 0) and (len(sp1_sp2_str) == 0):
            print 'neither'
            counts['neither_win'] += 1

        print >> hits_out, "%s:%s:%s-%s %s %s" % (sp1,sp1_chr,sp1_winstart,sp1_winend,sp1_sp1_str,sp1_sp2_str)
        print >> alns_out, sp1_sp2_aln

    print >> counts_out, "sp1_matrix: %s" % mtx1
    print >> counts_out, "sp2_matrix: %s" % mtx2
    print >> counts_out, "sp1: %s" % sp1
    print >> counts_out, "sp2: %s" % sp2
    print >> counts_out, "window: %s" % winsize
    for (k,v) in counts.items():
        print >> counts_out, "%s : %s" % (k,v)
        

def comp_site_aln(sp1_sites,sp1,sp2_sites,sp2,aln):
    sp1_sp2_count = 0
    for s1site in sp1_sites:
        for s2site in sp2_sites:
            print (s1site,s2site)
            s2start = aln.transformSpCoord(sp1,sp2,s1site['start'])
            s2end = aln.transformSpCoord(sp1,sp2,s1site['end'])
            if ((s2site['start'] < s2start + 3) and (s2site['start'] > s2start - 3) or
                (s2site['end'] < s2end + 3) and (s2site['end'] > s2end - 3)):
                sp1_sp2_count += 1
    sp1_only = len(sp1_sites) - sp1_sp2_count
    sp2_only = len(sp2_sites) - sp1_sp2_count
    print (sp1_sp2_count,sp1_only,sp2_only)
    return (sp1_sp2_count,sp1_only,sp2_only)

if __name__ == "__main__":
    main(sys.argv[1:])

