#!/usr/bin/env python

import sys
import os
import re
import numpy as np
import matplotlib.pyplot as plt

def main(args):
    
    if len(args) < 2:
        print "plot_domain_rates_from_files.py <rate_dir> <tf_name>"
        sys.exit(1)

    dir = args[0]
    dirs = os.listdir(dir)

    domains = {}
    target_domains = []
    target_fam = []
    medians = {}

    plot_family = False
    domain_nms = []

    species = None

    print dirs

    for d in dirs:
        dname_re = re.search("(\S+)_dom_rates",d)
        if not dname_re:
            continue
        dname = dname_re.group(1)
        fname = dname_re.group(0)
        rates = {}
        ratefile = open("%s%s/%s_pw_rates_domains" % (dir,fname,dname))
        species_ln = ratefile.readline()
        species = species_ln.split('\t')[1:-1]
        for l in ratefile:
            spl = l.split('\t')[:-1]
            pname = spl[0]
            for (i,x) in enumerate(spl[1:]):
                if (x == '-'):
                    continue
                #print rates.keys()
                #print species
                print species[i]
                print rates.keys()
                if species[i] in rates.keys():
                    print 'yes'
                    rates[species[i]].append(float(x))
                else:
                    rates[species[i]] = [float(x)]
            name_search = re.search(args[1],pname)
            if name_search:
                plot_family = True
                d_rates = {}
                (d_rates['name'],rs) = (pname,spl[1:])
                for (i,x) in enumerate(species):
                    #print (i,x,rs[i])
                    r = rs[i]
                    if r == '-':
                        continue
                    d_rates[x] = float(r)
                target_domains.append(d_rates)
            #print rates
        if plot_family:
            print rates.items()
            famplot(dname,rates,target_domains)
            plot_family = False
        for (sp,vals) in rates.items():
            avals = np.array(vals)
            sp_med = np.median(avals)
            if sp in medians.keys():
                print len(medians[sp])
                #print 'app ' + sp
                #print medians[sp]
                medians[sp].append(sp_med)
            else:
                print 'new ' + sp
                medians[sp] = [sp_med]
        domain_nms.append(dname)

    #print medians    

    medlist = []

    print domain_nms

    for (i,n) in enumerate(domain_nms):
        m = {'name': n}
        meds = []
        for s in species:
            if len(medians[s]) <= i:
                meds.append(-1)
            else:
                meds.append(medians[s][i])
        m['meds'] = meds
        medlist.append(m)

    medlist.sort(key=(lambda x: sum(x['meds']) ))
    rout = open("FlyTF_domain_rate_table.txt","w")

    print medlist

    for m in medlist:
        print >> rout, "%s&%s\\" % (m['name'],"&".join(map(str,m['meds'])))
        print "%s&%s\\" % (m['name'],"&".join(map(str,m['meds'])))

    for sp,meds in medians.items():
        plt.figure(1,figsize=(7,7))
        plt.hist(medians[sp],bins=np.arange(0,.5,.01))
        plt.xlabel("Per-family Median Pairwise %s-%s AA rate" % ('dmel',sp))
        plt.ylabel("Family Count")
        plt.ylim((0,40))
        plt.xlim((0,.5))
        #plt.show()
        plt.savefig("%s-dmel-%s_all-domain-median-rates.svg" % ('FlyTFs',sp),format="svg")
        plt.close()
        rout = open("%s-%s_domains_median-rates.txt" % ('dmel',sp),"w")
        print >> rout, "\t".join([str(x) for x in meds])


def famplot(name,rates,target_doms):
    for (species,rates) in rates.items():
        ratea = np.array(rates)
        med = np.median(ratea)
        print "%s %s Median: %f" % (species,name,med)
       
        fig = plt.figure(1,figsize=(7,7))
        p = plt.hist(ratea,bins=50)
    #    plt.semilogx()
        plt.xlabel("Pairwise %s-%s AA rate for %s" % ('dmel',species,name))
        plt.ylabel("Protein Count")
        
        print target_doms
        
        for (i,h) in enumerate(target_doms):
            if species not in h.keys():
                continue
            r = h[species]
            n = h['name']
            print h
            plt.annotate(n, xy=(r,1),  xycoords='data',
                         xytext=(50, 50*i+30), textcoords='offset points',
                         arrowprops=dict(arrowstyle="->")
                         )
        #plt.show()
        plt.savefig("%s_%s-%s_median-rates.svg" % (name,'dmel',species),format="svg")
        plt.close()

if __name__ == "__main__":
    main(sys.argv[1:])

 
    
                                   
