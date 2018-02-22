#!/usr/bin/env python
"""
patser_annotate_genome.py

search all chromosomes in a genome sequence file for specified matrix.  Annotate hits as rows in sqlite table
"""
import sys
import sqlite3
import fasta_subseq_2
import patser_tools
#from multiprocessing import Pool
#from pprint import pprint as pp

class searchObj(object):
    def __init__(self,
                 chrObj=None,
                 seq_name = None,
                 matrix=None,
                 matrix_name=None,
                 annotation=None):
        self.chrObj = chrObj
        self.seq = None
        self.seq_name = seq_name
        if chrObj:
            self.seq = chrObj['sequence'][0:(chrObj['sequence'].length - 1)]
        self.matrix = matrix
        self.matrix_name = matrix_name
        self.annotation = annotation
        print "created %s object" % chrObj['ID']
 
    def patSearch(self):
        annot=None
        try:
            annot = patser_tools.makePatserAnnotation(sequence=self.seq,matrix=self.matrix)
        except Exception:
           print "warning: Exception for seq %s" % (self.seq)
           annot = None
        self.annotation = annot

def search(s):
    print "starting search %s..." % s.chrObj["ID"]
    s.patSearch()
    if s.annotation:
        print "search complete: %s" % s.chrObj["ID"]
        return s
    else:
        print "search failed: %s!" % s.chrObj["ID"]
        return None

def _main(args):
    
    if len(args) != 4:
        print "usage: patser_annotate_genome_noxgrid.py <machi_db> <genome_seq> <matrix_file> <matrix_name>"
        sys.exit(0)

    #processes = int(args[4])

    # setup database
    conn = sqlite3.connect(args[0])
    conn.row_factory = sqlite3.Row
    cur = conn.cursor()

    cur.execute("""CREATE TABLE IF NOT EXISTS matrix (matrix_key INTEGER PRIMARY KEY,
                                                      name TEXT,
                                                      file TEXT)""")

    cur.execute("""CREATE TABLE IF NOT EXISTS patser_hit (patser_hit_key INTEGER PRIMARY KEY,
                                                          chr TEXT,
                                                          start INT,
                                                          end INT,
                                                          strand INT,
                                                          score FLOAT,
                                                          pval FLOAT,
                                                          matrix_key INT)""")

    cur.execute("""SELECT * FROM matrix WHERE file = ?""",(args[2],))
    matrix_exists = cur.fetchall()
    mtx_id = None
    if not matrix_exists:
        cur.execute("""INSERT INTO matrix VALUES (NULL,?,?)""", (args[3],args[2]))
        mtx_id = cur.lastrowid
    else:
        mtx_id = matrix_exists[0]["matrix_key"]


    # open fasta
    fasta = fasta_subseq_2.FastaDB()
    fasta.openFastaFile(args[1])
    
    jobs = []
    
    for (name,chr) in fasta.items():
        srch = searchObj(chrObj=chr,
                         seq_name = name,
                         matrix=args[2],
                         matrix_name=args[3])
        print srch
        jobs.append(srch)
    print jobs

    #pool = Pool(processes)
    
    #results = pool.imap(search,jobs)
    for j in jobs:
        s = search(j)
        print "inserting %s, %i tags" % (s.seq_name, len(s.annotation.getAllFeatures()))
        for feature in s.annotation.getAllFeatures():
            print >> sys.stderr, feature
            cur.execute("INSERT INTO patser_hit VALUES (NULL,?,?,?,?,?,?,?)",
                        (s.seq_name,feature.start,feature.end,feature.tags["strand"],feature.tags["score"],feature.tags["pval"],mtx_id))

    conn.commit()
    conn.close()
if __name__ == "__main__":
    _main(sys.argv[1:])
