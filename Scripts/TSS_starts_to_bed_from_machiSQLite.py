#!/usr/bin/env python
"""
get_TSS_starts_from_machiSQLite.py
"""

import sys
import sqlite3

def _main(args):
    
    if args != 1:
        print "usage: TSS_starts_to_bed_from_machiSQLite.py <sqlite3_DB> <window>"
        sys.exit(0)
        
    conn = sqlite3.connect(args[0])
    conn.row_factory = sqlite3.Row
    cur = conn.cursor()

    machi_out = open(args[0] + "_machiTSS" + args[1] + ".bed","w")
    fb_out = open(args[0] + "_fbTSS" + args[1] + ".bed","w")
    
    cur.execute("""SELECT * FROM machi_tag JOIN clusters ON machi_tag.machi_tag_key = clusters.peak_machi_tag_key""")
    machi_tss = cur.fetchall()
    cur.execute("""SELECT * from fb_tss""")
    fb_tss = cur.fetchall()
    
    for machi in machi_tss:
        if machi["strand"] == "-":
            print >> machi_out "%s\t%i\t%i\t%s\t%i\t%s" % (machi["chr"],(machi["end"] - int(args[1])),(machi["end"] + int(args[1])),machi["machi_tag_key"],machi["machi_tag.tags"],machi["strand"])
        else:
            print >> machi_out "%s\t%i\t%i\t%s\t%i\t%s" % (machi["chr"],(machi["start"] - int(args[1])),(machi["start"] + int(args[1])),machi["machi_tag_key"],machi["machi_tag.tags"],machi["strand"])

    for fb in fb_tss:
        if fb["strand"] == "-":
            print >> tss_out "%s\t%i\t%i\t%s\t.\t%s" % (fb["chr"],(fb["end"] - int(args[1])),(fb["end"] + int(args[1])),fb["fb_tss_key"],fb["strand"])
        else:
            print >> tss_out "%s\t%i\t%i\t%s\t.\t%s" % (fb["chr"],(fb["start"] - int(args[1])),(fb["start"] + int(args[1])),fb["fb_tss_key"],fb["strand"])
