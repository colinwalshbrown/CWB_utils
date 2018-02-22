#!/usr/bin/env python
"""
build table of orthologs from two flybase gffs
"""

import sys
import re
import sqlite3

def main(args):
    
    if len(args) < 2:
        print "usage: make-ortholog-table.py <GFF_SP1_REFERENCE> <NAME_SP1> <GFF2> <NAME_SP2_ORTHO> <sqlitedb_name>"
        sys.exit(1)

    print args
    sp1_name = args[1]
    sp2_name = args[3]
    sp1_gff = args[0]
    sp2_gff = args[2]

    conn = sqlite3.connect(args[4])
    conn.row_factory = sqlite3.Row
    cur = conn.cursor()

    cur.execute("DROP TABLE IF EXISTS %s_gene" % sp1_name)
    cur.execute("DROP TABLE IF EXISTS %s_gene" % sp2_name)
    cur.execute("DROP TABLE IF EXISTS %s_%s_orth" % (sp1_name,sp2_name))

    cur.execute("""CREATE TABLE IF NOT EXISTS %s_gene (%s_gene_key INTEGER PRIMARY KEY AUTOINCREMENT,
                                                     chr TEXT,
                                                     start INT,
                                                     end INT,
                                                     strand TEXT,
                                                     FBid TEXT,
                                                     name TEXT,
                                                     ann_str TEXT)""" % (sp1_name,sp1_name))
    cur.execute("""CREATE INDEX %s_gene_idx ON %s_gene (chr,start,end)""" % (sp1_name,sp1_name))

    cur.execute("""CREATE TABLE IF NOT EXISTS %s_gene (%s_gene_key INTEGER PRIMARY KEY AUTOINCREMENT,
                                                     chr TEXT,
                                                     start INT,
                                                     end INT,
                                                     strand TEXT,
                                                     FBid TEXT,
                                                     name TEXT,
                                                     ann_str TEXT)""" % (sp2_name,sp2_name))
    cur.execute("""CREATE INDEX %s_gene_idx ON %s_gene (chr,start,end)""" % (sp2_name,sp2_name))

    
    cur.execute("""CREATE TABLE IF NOT EXISTS %s_%s_orth (%s_gene_key INT,
                                                         %s_gene_key INT)""" % (sp1_name,sp2_name,sp1_name,sp2_name))

    for gff_line in open(sp1_gff):
        # cycle GFF lines
        gff_split = gff_line[:-1].split("\t")
        if (re.match("^#",gff_split[0]) or (len(gff_split) < 8)):
            if (re.match("^##FASTA",gff_split[0])):
                break
            else:
                continue
        
        if not gff_split[2] == 'gene':
            continue
        ID_match = re.search("ID=(\S+?);",gff_split[8])
        Name_match = re.search("Name=(\S+?);",gff_split[8])
        cur.execute("INSERT INTO %s_gene VALUES (?,?,?,?,?,?,?,?)" % sp1_name,(None,gff_split[0],gff_split[3],gff_split[4],gff_split[6],ID_match.group(1),Name_match.group(1),gff_split[8]))
        
    print "Done reading %s" % sp1_gff
    conn.commit()

    for gff_line in open(sp2_gff):
        # cycle GFF lines
        gff_split = gff_line[:-1].split("\t")

        if (re.match("^#",gff_split[0]) or (len(gff_split) < 8)):
            if (re.match("^##FASTA",gff_split[0])):
                break
            else:
                continue

        if not gff_split[2] == 'gene':
            continue

        ID_match = re.search("ID=(\S+?);",gff_split[8])
        Name_match = re.search("Name=(\S+?);",gff_split[8])
        cur.execute("INSERT INTO %s_gene VALUES (?,?,?,?,?,?,?,?)" % sp2_name,(None,gff_split[0],gff_split[3],gff_split[4],gff_split[6],ID_match.group(1),Name_match.group(1),gff_split[8])) 

    print "Done reading %s" % sp2_gff
    conn.commit()

    count = 0

    for gff_line in open(sp2_gff):
        # cycle GFF lines
        gff_split = gff_line[:-1].split("\t")

        if (re.match("^#",gff_split[0]) or (len(gff_split) < 8)):
            if (re.match("^##FASTA",gff_split[0])):
                break
            else:
                continue

        if ((not gff_split[2] == 'orthologous_to') or
            (not gff_split[1] == 'FlyBase')):
            continue
        
        ID_match = re.search("ID=(\S+?)_\S+?;",gff_split[8])
        to_species_match = re.search("to_species=(\S+?);",gff_split[8])
        to_id_match = re.search("to_id=(\S+?);",gff_split[8])

#        print gff_split[8]
#        print ID_match.group(1)
#        print to_species_match.group(1)
#        print to_id_match.group(1)

        if not to_species_match.group(1).lower() == sp1_name.lower():
            continue
        
        cur.execute("SELECT * from %s_gene WHERE FBid = ?" % (sp2_name,),(ID_match.group(1),))
        ID_result = cur.fetchone()
        cur.execute("SELECT * from %s_gene WHERE FBid = ?" % (sp1_name,),(to_id_match.group(1),))
        to_ID_result = cur.fetchone()

        if not to_ID_result:
            print gff_split[8]
            targ_coord_match = re.search("Target=(.+)",gff_split[8])
            if not targ_coord_match:
                print >> sys.stderr, "WARNING: no %s_gene entry found for %s" % (sp1_name,gff_split[8])
                continue
                
            print targ_coord_match.group(1).split()[0:3]
            (chr,start,end) = targ_coord_match.group(1).split()[0:3]
            print "%s %s %s" % (chr,start,end)
            cur.execute("SELECT %s_gene_key from %s_gene where (chr = ? AND (start = ? AND end = ?))" % (sp1_name,sp1_name),(chr,start,end))
            to_ID_result = cur.fetchone()
            if not to_ID_result:
                print >> sys.stderr, "WARNING: no %s_gene entry found for %s" % (sp1_name,gff_split[8])
                continue
            else:
                print >> sys.stderr, "%s mapped to %s" % (gff_split[8],to_ID_result)

        cur.execute("INSERT INTO %s_%s_orth VALUES (?,?)" % (sp1_name,sp2_name),(to_ID_result[sp1_name+"_gene_key"],ID_result[sp2_name+"_gene_key"]))
        count += 1

    print "Done building orthology table; %d orthologs inserted" % count

    conn.commit()
    conn.close()

if __name__ == "__main__":
    main(sys.argv[1:])
