#!/usr/bin/env python

import sys
import re
import patser_tools
import math
import fasta_subseq_2

def _main(args):

    if len(args) < 1:
        print "usage: patser_list.py [mtx1][mtx2] ... < seqs.fa"
        sys.exit(1)
    
    matrices = args
    seqs = sys.stdin
    hits = {}
    name = ""
    seq = ""
    pssms = convertFreqMtx(matrices)
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
        print "%s_score\t%s_pval\t%s_PSSM_score" % (x,x,x),
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
            pssm_scores = [scoreSeq(pssms[x],h['seq']),scoreSeq(pssms[x],fasta_subseq_2.revcomp(h['seq']))]
            score = None
            if pssm_scores[1] > pssm_scores[0]:
                score = pssm_scores[1]
            else:
                score = pssm_scores[0]
            print str(score) + "\t",
        print ""

def scoreSeq(mtx,seq):
    score = 0
    for (x,fn) in enumerate(mtx):
        try:
            pos_scr = fn[seq[x].upper()]
            score += pos_scr
        except KeyError:
            print "Error: lengths of matrix (%d) and seq %s (%d) don't match..." % (len(mtx),seq,len(seq))
            return None
    return score

def convertFreqMtx(matrix_files):
    # assumes patser vertical matrix format
    pssms = {}
    for m in matrix_files:
        mtx_name = ".".join(m.split("/")[-1].split(".")[:-1])
        mtx = []
        mtx_order = None
        for m_line in open(m):
            #print m_line
            mtx_re = re.search("\s*(\S+)\s*\|\s*(\S+)\s*\|\s*(\S+)\s*\|\s*(\S+)\s*\|.*",m_line)
            mtx_counts = None
            if mtx_re and mtx_order:
                mtx_counts = map(int,mtx_re.group(1,2,3,4))           
            elif mtx_re and not mtx_order:
                mtx_order = mtx_re.group(1,2,3,4)
                continue
            else:
                continue
            mtx_probs = {}
            bg_prob = {'A':.25, # could fix for completeness...
                       'T':.25,
                       'G':.25,
                       'C':.25}
            nsites = sum(mtx_counts)
            ps_cnt = math.sqrt(nsites) # could be something else...see Wasserman(2004) for refs
            mtx_fn = lambda ns,ps,cn,nt: math.log((((cn + ps) / (ns + 4*ps))/bg_prob[nt]),2)
            for (n,x) in enumerate(mtx_order):
                score = mtx_fn(nsites,ps_cnt,mtx_counts[n],x)
                mtx_probs[x] = score
            mtx.append(mtx_probs)
        pssms[mtx_name] = mtx
    #print pssms
    return pssms

if __name__ == "__main__":
    _main(sys.argv[1:])
