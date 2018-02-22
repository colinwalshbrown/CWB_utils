#!/usr/bin/env python

import sys
import re
import patser_tools

def _main(args):

    if len(args) < 1:
        print "usage: patser_list.py [mtx1][mtx2] ... < seqs.fa"
        sys.exit(1)
    
    matrices = args
    seqs = sys.stdin
    hits = {}
    name = ""
    seq = ""
    for s in seqs:
        nameres = re.search(">(\S+)",s)

        if nameres and not (name == ""):
            hits[name] = {'seq' : seq}
            name = nameres.group(1)
            seq = ""
        elif nameres and (name == ""):
            #print "1"
            name = nameres.group(1)
            seq = ""
        else:
            seq += s[:-1]
    if not (name == ""):
        hits[name] = {'seq':seq}

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
                print >> sys.stderr, "Sequence %s: No hit for matrix %s in %s" % (name,mtx,d['seq'])
                continue
            #print hit
            d[hit.tags['motif_name']] = hit
            if hit.tags['motif_name'] not in mtx_names:
                mtx_names.append(hit.tags['motif_name'])
    #print hits
    print "name\tsequence\t",
    for x in mtx_names:
        print "%s_score\t%s_pval\t" % (x,x),
    print ""
    for (name,h) in hits.iteritems():
        matrices = [x for x in h.keys() if not x == 'seq']
        print "%s\t%s\t" % (name,h['seq']) ,
        for x in matrices:
            print str(h[x].tags['score']) + "\t",
            if 'pval' in h[x].tags.keys():
                print str(h[x].tags['pval']) + "\t",
            else:
                print "-",
        print ""

if __name__ == "__main__":
    _main(sys.argv[1:])
