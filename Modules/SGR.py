#!/usr/bin/env python

import sys
import os
import re
import pickle
import sqlite3
import alignment

DB_DIR = os.environ["HOME"] + "/scripts/python_modules/SGR_DBs/"

class SGR():
    
    def __init__(self,sgrfile,reindex=False):
        self.file = sgrfile
        self.sgrDBs = None
        self.DBconn = None
        self.DBcur = None
        self.DBname = None
        try:
            self.sgrDBs = pickle.load(open(os.environ["HOME"] + "/scripts/python_modules/SGR_DB_INDEXES"))
        except:
            print >> sys.stderr, "Building new SGR index file..."
            self.sgrDBs = {}    
        if (not self.sgrDBs) or (sgrfile not in self.sgrDBs.keys()):
            print >> sys.stderr, "%s not in index file, building database..." % (sgrfile,)
            self._ParseSGR(sgrfile)
            self.sgrDBs[sgrfile] = self.DBname
            self.sgrDBs = pickle.dump(self.sgrDBs,open(os.environ["HOME"] + "/scripts/python_modules/SGR_DB_INDEXES","w"))
        elif reindex:        
            print >> sys.stderr, "Reindexing %s..." % (sgrfile,)
            os.remove(self.sgrDBs[sgrfile])
            self._ParseSGR(sgrfile)
            self.sgrDBs[sgrfile] = self.DBname
            self.sgrDBs = pickle.dump(self.sgrDBs,open(os.environ["HOME"] + "/scripts/python_modules/SGR_DB_INDEXES","w"))
        elif sgrfile in self.sgrDBs.keys():
            self.DBconn = sqlite3.connect(self.sgrDBs[sgrfile])
            self.DBcur = self.DBconn.cursor()
            self.DBname = self.sgrDBs[sgrfile]
            self.DBcur.row_factory = sqlite3.Row
            

    def _ParseSGR(self,file):
        DB_name = ".".join(file.split("/")[-1]) + "_SGR_DB"
        self.DBconn = sqlite3.connect(DB_DIR + DB_name)
        self.DBconn.row_factory = sqlite3.Row
        self.DBcur = self.DBconn.cursor()
        self.DBname = DB_DIR + DB_name
        self.DBcur.execute("""DROP TABLE IF EXISTS sgr_values""")
        self.DBcur.execute("""CREATE TABLE IF NOT EXISTS sgr_values (sgr_values_key INTEGER PRIMARY KEY AUTOINCREMENT,
                                                                     chr TEXT,
                                                                     loc INT,
                                                                     val FLOAT)""")
        self.DBcur.execute("""CREATE INDEX IF NOT EXISTS sgr_values_idx ON sgr_values (chr,loc)""")
        for l in open(file):
            sgr_re = re.search("(\S+)\s+(\S+)\s+(\S+)",l)
            chr = sgr_re.group(1)
            loc = int(sgr_re.group(2))
            val = float(sgr_re.group(3))
            
            self.DBcur.execute("""INSERT INTO sgr_values VALUES (?,?,?,?)""",(None,chr,loc,val))
        self.DBconn.commit()

    def SGRslice(self,chr,start,end):
        val_dict = {}
        self.DBcur.execute("""SELECT * FROM sgr_values WHERE ((chr = ?) AND (loc BETWEEN ? and ?))""",(chr,start,end))
        slice = self.DBcur.fetchall()
        val_dict['step'] = slice[1]['loc'] - slice[0]['loc']
        val_dict['chr'] = chr
        val_dict['start'] = min(map(lambda x:x['loc'],slice))
        val_dict['end'] = max(map(lambda x:x['loc'],slice))
        sorted_slice = sorted(slice,key=(lambda x:x['loc']))
        val_dict['vals'] = map(lambda x:x['val'],sorted_slice)
        return val_dict
                                                                 
