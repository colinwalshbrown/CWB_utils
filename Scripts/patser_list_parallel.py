#!/usr/bin/env python

import sys
import re
import patser_tools
import xgrid_tools

CHUNKSIZE=10000

class seq_scanner():
    def __init__(self,
                 seqs=None,
                 matrices=None):
        self.seqs=seqs
        self.matrices=matrices
        self.results=None
        self.header=None
    
    def runXgrid(self):
        (self.results,self.header) = scan_seqs(self.seqs,self.matrices)
        return(self.results,self.header)

def _main(args):
    
    parallel = False
    if len(args) < 1:
        print "usage: patser_list.py --parallel [mtx1][mtx2] ... < seqs.fa"
        sys.exit(1)
    if "-p" in args:
        parallel = True
        args.remove("-p")
    elif "--parallel" in args:
        parallel = True
        args.remove("--parallel")

    matrices = args
    seqs = sys.stdin
    if not parallel:
        (header,result) = scan_seqs(seqs,matrices)
        print header
        print results
    else:
        seqlist = []
        jobs 
        for (i,s) in enumerate(seqs):
            if (i % CHUNKSIZE) ==                                   
    
    

def scan_seqs(seqs,matrices):
    hits = {}
    name = ""
    seq = ""
    header = ""
    results = ""
    for s in seqs:
        nameres = re.search(">(\S+)",s)

        if nameres and not (name == ""):
            hits[name] = {'seq' : seq}
            name = nameres.group(0)
            seq = ""
        elif nameres and (name == ""):
            #print "1"
            name = nameres.group(0)
            seq = ""
        else:
            seq += s[:-1]

    #print hits
    mtx_names = []
    for (name,d) in hits.iteritems():
        for mtx in matrices:
            hit_annot = patser_tools.makePatserAnnotation(sequence=d['seq'],matrix=mtx,seqname=name,scorecut=-100)
            features = hit_annot.getAllFeatures()
            hit = None
            if len(features) > 0:
                max = features[0]
                for x in features:
                    if x.tags['score'] > max.tags['score']:
                        max = x
                hit = max
            else:
                print >> sys.stderr, "No hit for matrix %s in %s" % (mtx,d['seq'])
                continue
            #print hit
            d[hit.tags['motif_name']] = hit
            if hit.tags['motif_name'] not in mtx_names:
                mtx_names.append(hit.tags['motif_name'])
    #print hits
    header = "name\tsequence\t",
    for x in mtx_names:
        header += "%s_score\t%s_pval\t" % (x,x),
    for (name,h) in hits.iteritems():
        matrices = [x for x in h.keys() if not x == 'seq']
        result += "%s\t%s\t" % (name,h['seq']) ,
        for x in matrices:
            result += str(h[x].tags['score']) + "\t",
            if 'pval' in h[x].tags.keys():
                result += str(h[x].tags['pval']) + "\t",
            else:
                result += "-",
    return (header,result)

if __name__ == "__main__":
    _main(sys.argv[1:])
