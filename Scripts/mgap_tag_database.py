#!/usr/bin/env python

import sys
import sqlite3
import re

def _main(args):
    
    if len(args) != 3:
        print "usage: mgap_tag_database.py <sqlite> <mgap_map_file> <genome_name>"
        sys.exit(1)

    con = sqlite3.Connection(args[0])
    cur = con.cursor()

    cur.execute("""CREATE TABLE IF NOT EXISTS %s_mgap_tags (%s_mgap_tags_key INTEGER PRIMARY KEY AUTOINCREMENT,
                                                       chr TEXT,
                                                       start INT,
                                                       read_strand TEXT,
                                                       called_strand TEXT,
                                                       chr_strand TEXT,
                                                       percent_id FLOAT,
                                                       mismatches INT,
                                                       indels INT,
                                                       template TEXT)""" % (args[2],args[2]))
    cur.execute("""CREATE INDEX IF NOT EXISTS %s_mgap_tags_idx ON %s_mgap_tags (chr,start)""" % (args[2],args[2]))
    curtag = None
    for i in open(args[1]):
        head_srch = re.search("^>.+TEMPLATE\:\s+(\S+)\s+DIRECTION\:\s+(\S+)",i)
        pth_srch = re.search("Paths\s+\((\S+)\)\:",i)
        dir_srch = re.search("cDNA\s+direction\:\s+(\S+)",i)
        loc_srch = re.search("Path.+?genome\s+(\S+)\:(\S+?)\.\.\S+\s+\((\S+)\s+bp\)",i)
        pid_srch = re.search("Percent identity\:\s+(\S+)\s+.+?matches\,\s+(\S+)\s+mismatches\,\s+(\S+)\s+indels",i)
        if (head_srch) and (not curtag):
            curtag = {'template' : head_srch.group(1),
                      'read_strand' : head_srch.group(2)}
        elif (not curtag):
            continue
        elif (head_srch):
            print curtag
            cur.execute("""INSERT INTO %s_mgap_tags VALUES (?,?,?,?,?,?,?,?,?,?)""" % args[2],(None,curtag['chr'],curtag['start'],curtag['read_strand'],
                                                                                     curtag['called_strand'],curtag['chr_strand'],curtag['percent_id'],
                                                                                     curtag['mismatches'],curtag['indels'],curtag['template']))
            curtag = {'template' : head_srch.group(1),
                      'read_strand' : head_srch.group(2)}
        elif (pth_srch):
            if pth_srch.group(1) != '1':
                curtag = None
        elif (dir_srch):
            curtag['called_strand'] = dir_srch.group(1)
        elif (loc_srch):
            curtag['chr'] = loc_srch.group(1)
            curtag['start'] = "".join(loc_srch.group(2).split(","))
            if int(loc_srch.group(3)) < 0:
                curtag['chr_strand'] = "-"
            else:
                curtag['chr_strand'] = "+"
        elif (pid_srch):
            curtag['percent_id'] = pid_srch.group(1)
            curtag['mismatches'] = pid_srch.group(2)
            curtag['indels'] = pid_srch.group(3)
    
    con.commit()
    con.close()
        

if __name__ == "__main__":
    _main(sys.argv[1:])
