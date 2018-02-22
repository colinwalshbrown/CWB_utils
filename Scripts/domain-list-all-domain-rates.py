#!/usr/bin/env python

import sys
import sqlite3
import matplotlib.pyplot as plt
import numpy as np

def main(args):
    
    if len(args) < 2:
        print "usage: domain-db-all-domain-rates.py <db> <species> <output_name> < <domain_list>"
        sys.exit(1)

    conn = sqlite3.connect(args[0])
    conn.row_factory = sqlite3.Row
    cur = conn.cursor()
    sp1 = 'dmel'
    sp2 = args[1]
    out = args[2]

    domains = {}

    dom_names = sys.stdin
    dom_list = []

    for n in dom_names:
        cur.execute("""SELECT * FROM pfam_model WHERE name = ?""", (n[:-1],))
        dom = cur.fetchone()
        if not dom:
            print "no %s domain model!" % (n,)
            continue
        d = {'desc':dom['description'],
             'length':dom['length'],
             'acc':dom['accession'],
             'pfam_id':dom['pfam_id']}
        domains[n] = d

    print "%s domains found" % len(domains.keys())
    
    
    for (k,d) in domains.items():
        cur.execute("""SELECT r.rate 
                       FROM domain_pw_aa_rates r JOIN dmel_pfam_domain d 
                       ON r.dmel_pfam_domain_id = d.dmel_pfam_domain_id 
                       WHERE ((d.pfam_id = ?) AND (r.species1 = ?)) AND (r.species2 = ?)""", (str(d['pfam_id']),sp1,sp2))
        rates_result = cur.fetchall()
        rates = np.array([float(x['rate']) for x in rates_result])
        print rates
        med = np.median(rates)
        d['med'] = med
        print "%s: %s median = %s" % (k,'dmel_'+sp2,med)
    
    median_rates = [x['med'] for x in domains.values() if x['med'] >= 0]
    plt.figure(1,figsize=(7,7))
    plt.hist(median_rates,bins=50)
    plt.xlabel("Per-family Median Pairwise %s-%s AA rate" % (sp1,sp2))
    plt.ylabel("Family Count")
    plt.show()
    plt.savefig("%s_%s-%s_domain-median-rates.svg" % (out,sp1,sp2),format="svg")
    
    rout = open(out + "_median-rates.txt","w")
    print >> rout, "\t".join([str(x) for x in median_rates])
                 
if __name__ == "__main__":
    main(sys.argv[1:])
