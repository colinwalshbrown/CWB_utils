#!/usr/bin/env python

import sys
import sqlite3
import matplotlib.pyplot as plt
import numpy as np

def main(args):
    
    if len(args) < 2:
        print "usage: domain-db-all-domain-rates.py <db> <species>"
        sys.exit(1)

    conn = sqlite3.connect(args[0])
    conn.row_factory = sqlite3.Row
    cur = conn.cursor()
    sp1 = 'dmel'
    sp2 = args[1]

    domains = {}

    cur.execute("""SELECT * FROM pfam_model""")
    dom_list = cur.fetchall()
    
    for x in dom_list:
        d = {'desc':x['description'],
             'length':x['length'],
             'acc':x['accession'],
             'pfam_id':x['pfam_id']}
        domains[x['name']] = d

    print "%s domains found" % len(domains.keys())
    
    cnt = 0
    for (k,d) in domains.items():
        print "searching %s %s..." % (k,d['pfam_id'])
        cur.execute("""SELECT r.rate 
                       FROM domain_pw_aa_rates r JOIN dmel_pfam_domain d 
                       ON r.dmel_pfam_domain_id = d.dmel_pfam_domain_id 
                       WHERE ((d.pfam_id = ?) AND (r.species1 = ?)) AND (r.species2 = ?)""", (str(d['pfam_id']),sp1,sp2))
        rates_result = cur.fetchall()
        print len(rates_result)
        rates = np.array([float(x['rate']) for x in rates_result])
        med = np.median(rates)
        d['med'] = med
        print "%d %s: %s median = %s" % (cnt,k,'dmel_' + sp2,med)
        cnt += 1
    
    median_rates = np.array([x['med'] for x in domains.values() if x['med'] >= 0])
    print median_rates
    plt.figure(1,figsize=(7,7))
    plt.hist(median_rates,bins=50)
    plt.xlabel("Per-family Median Pairwise %s-%s AA rate" % (sp1,sp2))
    plt.ylabel("Family Count")
    plt.show()
    plt.savefig("%s_all-domain-median-rates.svg" % (args[0],),format="svg")
    rout = open("%s-%s_domains_median-rates.txt" % (sp1,sp2),"w")
    print >> rout, "\t".join([str(x) for x in median_rates])
                 
if __name__ == "__main__":
    main(sys.argv[1:])
