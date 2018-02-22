#!/usr/bin/env python

import sys
import sqlite3
import fasta_subseq_2
import patser_tools
from pprint import pprint as pp

class searchObj():
    def __init__(self,
                 seq_obj=None,
                 matrix=None,
                 seq_name=None,
                 matrix_name=None,
                 annotation=None):
        self.seq_obj = seq_obj
        self.seq_name = seq_name
        self.matrix = matrix
        self.matrix_name = matrix_name
        self.annotation = annotation
 
    def patSearch(self):
        annot=None
        #try:
        annot = patser_tools.makePatserAnnotation(sequence=self.seq,matrix=self.matrix)
        #except Exception, e:
        #    print "warning: %s for seq %s" % (type(e),self.seq_name)
        #    annot = None
        self.annotation = annot
        
    def runXgrid(self):
        self.seq = self.seq_obj[0:(self.seq_obj.length - 1)]
        self.patSearch()

class insertObj():
    def __init__ (self,
                machi_db=None,
                matrix_key=None,
                annotations=None):
        self.machi_db = machi_db
        
        if machi_db:
            self.machi_conn = sqlite3.connect(self.machi_db)
            self.machi_conn.row_factory = sqlite3.Row
            self.machi_cur = self.machi_conn.cursor()
        
        self.matrix_key = matrix_key
        self.annotations = annotations

    def runXgrid(self):
        
        for (chr,ann) in annotations.items():
            for feature in ann.getAllFeatures():
                self.machi_cur.execute("INSERT INTO patser_hit VALUES (NULL,?,?,?,?,?,?,?)",
                                       (chr,feature.start,feature.end,feature.tags["strand"],feature.tags["score"],feature.tags["pval"],matrix_key))

        self.machi_conn.commit()
        self.machi_conn.close()
