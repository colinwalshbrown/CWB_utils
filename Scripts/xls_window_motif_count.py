#!/usr/bin/env python

import sys
import sqlite3
import numpy as np
import matplotlib.pyplot as plt

BINSIZE = 10

def _main(args):
    
    if len(args) != 3:
        print ("usage: mgap_txn_start_motif_count.py <mgap_motif_db> <xls> <window>")
        sys.exit(1)

    con = sqlite3.connect(args[0])
    con.row_factory = sqlite3.Row
    cur = con.cursor()

    cur.execute("""select * from matrix""")
    matrices = cur.fetchall()

    win = int(args[2])
    hit_cnt = {}
    for m in matrices:
        hit_cnt[m['name']] = 0

    total = 0
    dir = {}
    for ln in open(args[1]):
        total += 1
        sp = ln[:-1].split()
        pk = int(sp[1]) + int(sp[4])
#         cur.execute("""SELECT * FROM patser_hit WHERE ((matrix_key = ?) AND ((chr = ? ) AND (start BETWEEN ? AND ?)))""",('dpse_her',sp[0],pk - win,pk + win))
#         dp_hits = cur.fetchall()
#         dir = 0
#         if len(dp_hits) > 1:
#             if dp_hits[0]['strand'] == '+':
#                 dir = 1
#             else:
#                 dir = -1
                
        for m in matrices:
            cur.execute("""SELECT * FROM patser_hit WHERE ((matrix_key = ?) AND ((chr = ?) AND (start BETWEEN ? AND ?)))""",(m['matrix_key'],sp[0],pk - win,pk + win))
            if len(cur.fetchall()) > 0:
                hit_cnt[m['name']] += 1
#            for hit in cur.fetchall():
                
#                 if dir > -1:
#                     hit_cnt[m['name']][(hit['start'] - (pk - win))/BINSIZE] += 1
#                 else:
#                     hit_cnt[m['name']][-(hit['start'] - (pk - win))/BINSIZE] += 1

    for (name,cnt) in hit_cnt.items():
        print "%s\t%d\t%f" % (name,cnt,float(cnt)/total)

         

if __name__ == "__main__":
    _main(sys.argv[1:])
