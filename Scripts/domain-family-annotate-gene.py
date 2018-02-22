#!/usr/bin/env python

import sys
import sqlite3
import matplotlib.pyplot as plt
import numpy as np

def main(args):
    
    if len(args) < 2:
        print "usage: domain-family-annotate-gene.py <db> <species> <domain_name> <protein_FBid>"
        sys.exit(1)

    conn = sqlite3.connect(args[0])
    conn.row_factory = sqlite3.Row
    cur = conn.cursor()
    sp1 = 'dmel'
    sp2 = args[1]
    fbID = args[3]
    dname = args[2]

    cur.execute("""SELECT * FROM pfam_model WHERE name = ?""", (dname,))
    dom = cur.fetchone()
    if len(dom) < 1:
        print "no %s domain model!" % (n,)
    d = {'desc':dom['description'],
         'length':dom['length'],
         'acc':dom['accession'],
         'pfam_id':dom['pfam_id']}
    domain = d
    
    cur.execute("""SELECT r.rate,g.FBpp,d.idx 
                   FROM domain_pw_aa_rates r JOIN dmel_pfam_domain d 
                   ON r.dmel_pfam_domain_id = d.dmel_pfam_domain_id JOIN dmel_gene g
                   ON d.dmel_gene_id = g.dmel_gene_id
                   WHERE ((d.pfam_id = ?) AND (r.species1 = ?)) AND (r.species2 = ?)""", (str(d['pfam_id']),sp1,sp2))
    rates_result = cur.fetchall()
    rates = np.array([float(x['rate']) for x in rates_result])
    med = np.median(rates)
    d['med'] = med
    print "%s %s: %s median = %s" % (d['pfam_id'],dname,sp1+"_"+sp2,med)

#    median_rates = [x['med'] for x in domains.values()]

    gene_hits = [x for x in rates_result if (x['FBpp'] == fbID)]
    print rates
    plt.figure(1,figsize=(7,7))
    p = plt.hist(rates,bins=50)
#    plt.semilogx()
    plt.xlabel("Pairwise %s-%s AA rate for %s" % (sp1,sp2,dname))
    plt.ylabel("Protein Count")

    print gene_hits

    for h in gene_hits:
        print h
        plt.annotate(h[1]+"-"+h[2], xy=(h['rate'],1),  xycoords='data',
                     xytext=(-50, 30), textcoords='offset points',
                     arrowprops=dict(arrowstyle="->")
                     )
    plt.show()
    plt.savefig("%s_%s-%s_median-rates.svg" % (args[0],sp1,sp2),format="svg")
    
                 
if __name__ == "__main__":
    main(sys.argv[1:])
