#!/usr/bin/env python

import sys
import sqlite3

def _main(args):
    
    if len(args) != 3:
        print ("usage: mgap_txn_start_motif_count.py <mgap_motif_db> <genome> <window>")
        sys.exit(1)

    con = sqlite3.connect(args[0])
    con.row_factory = sqlite3.Row
    cur = con.cursor()

    cur.execute("""SELECT * FROM matrix""")
    
    matrices = cur.fetchall()
    win = int(args[2])

    cur.execute("""SELECT * FROM %s_mgap_clusters""" % args[1])

    for cl in cur.fetchall():
        print "=" * 40
        print cl
        for m in matrices:
            cur.execute("""SELECT * FROM patser_hit WHERE ((matrix_key = ?) AND ((chr = ?) AND (start BETWEEN ? AND ?)))""",(m['matrix_key'],cl['chr'],(cl['peak'] - win),(cl['peak'] + win)))
            print "%s => " % m['name'],
            for hit in cur.fetchall():
                print (cl['peak'] - hit['start']),
                print " ",
            print ""

if __name__ == "__main__":
    _main(sys.argv[1:])
