#!/usr/bin/env python
"""
patser_annotate_genome.py

search all chromosomes in a genome sequence file for specified matrix.  Annotate hits as rows in sqlite table
"""
import sys
import sqlite3
import fasta_subseq_2
import patser_tools
import multiprocessing as mp
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
        print_annotation(s)
        return s
    else:
        print "search failed: %s!" % s.chrObj["ID"]
        return None

def print_annotation(sobj):
    outfile = open(sobj.seq_name + "_" + sobj.matrix_name + "_patsearch.bed","w")
    for feature in sobj.annotation.getAllFeatures():
        outln = "%s\t%d\t%d\t%s\t%f\t%s\t%f" % (sobj.seq_name,feature.start,feature.end,sobj.matrix_name + "_hit",feature.tags["score"],feature.tags["strand"],feature.tags["pval"])
        print >> outfile, outln
    outfile.close()

def _main(args):
    
    if len(args) != 3:
        print "usage: patser_annotate_genome_noxgrid.py <genome_seq> <matrix_file> <matrix_name>"
        sys.exit(0)

    # open fasta
    fasta = fasta_subseq_2.FastaDB()
    fasta.openFastaFile(args[0])
    
    jobs = []
    
    for (name,chr) in fasta.items():
        srch = searchObj(chrObj=chr,
                         seq_name = name,
                         matrix=args[1],
                         matrix_name=args[2])
        print srch
        jobs.append(srch)
    print jobs

    pool = mp.Pool()
    
    results = pool.map(search,jobs)

if __name__ == "__main__":
    _main(sys.argv[1:])
