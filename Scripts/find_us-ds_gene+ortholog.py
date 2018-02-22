#!/usr/bin/env python
"""

"""

import sys
import re
import sqlite3

def main(args):

    if len(args) < 3:
        print "usage: find_nearest_gene+ortholog.py <orth_db> <sp1_xls> <sp1_name> <sp2_xls> <sp2_name> <window>"
        sys.exit(1)

    conn = sqlite3.connect(args[0])
    conn.row_factory = sqlite3.Row
    cur = conn.cursor()

    xls1 = args[1]
    xls2 = args[3]
    sp1_name = args[2]
    sp2_name = args[4]
    win = int(args[5])

    xls1_d = []
    xls2_d = []
    
    for line in open(xls1):
        if re.match("#",line):
            continue
        spl = line[:-1].split()
        for pk in spl[4].split(',')[:-1]:
            entr = { 'chr' : spl[1],
                     'loc' : int(pk)}
            xls1_d.append(entr)
    
    for line in open(xls2):
        if re.match("#",line):
            continue
        spl = line[:-1].split()
        for pk in spl[4].split(',')[:-1]:
            entr = { 'chr' : spl[1],
                     'loc' : int(pk)}
            xls2_d.append(entr)

    for x in xls1_d:
        cur.execute("""select * from %s_gene where chr = ? and 
                                                 ((start between ? and ?) or
                                                 (end between ? and ?))""" % sp1_name, (x['chr'],x['loc']-win,x['loc']+win,x['loc']-win,x['loc']+win))
        hits = cur.fetchall()
        (us_hit,us_dist) = get_us_hit(x,hits,cur)
        (ds_hit,ds_dist) = get_ds_hit(x,hits,cur)

        if not us_hit:
            print >> sys.stderr, "No upstream gene found for %s" % x
            continue

        cur.execute("""SELECT %s_gene.* from %s_gene JOIN %s_%s_orth USING (%s_gene_key)
                                                     JOIN %s_gene USING (%s_gene_key)
                                                     WHERE %s_gene.%s_gene_key = ?""" % (sp2_name,sp1_name,sp1_name,sp2_name,sp1_name,sp2_name,sp2_name,sp1_name,sp1_name),(us_hit[sp1_name+"_gene_key"],))

        x['us_orth'] = cur.fetchone()
        x['us_gene'] = us_hit
        x['us_dist'] = -us_dist

        if not ds_hit:
            print >> sys.stderr, "No downstream gene found for %s" % x
            continue

        cur.execute("""SELECT %s_gene.* from %s_gene JOIN %s_%s_orth USING (%s_gene_key)
                                                     JOIN %s_gene USING (%s_gene_key)
                                                     WHERE %s_gene.%s_gene_key = ?""" % (sp2_name,sp1_name,sp1_name,sp2_name,sp1_name,sp2_name,sp2_name,sp1_name,sp1_name),(ds_hit[sp1_name+"_gene_key"],))

        x['ds_orth'] = cur.fetchone()
        x['ds_gene'] = ds_hit
        x['ds_dist'] = ds_dist
 
    for x in xls2_d:
        cur.execute("""select * from %s_gene where chr = ? and 
                                                 ((start between ? and ?) or
                                                 (end between ? and ?))""" % sp2_name, (x['chr'],x['loc']-win,x['loc']+win,x['loc']-win,x['loc']+win))
        hits = cur.fetchall()
        (us_hit,us_dist) = get_us_hit(x,hits,cur)
        (ds_hit,ds_dist) = get_ds_hit(x,hits,cur)
        #get_min_orth = get_min_orth(min_hit,cur)
        if not us_hit:
            print >> sys.stderr, "No upstream gene found for %s" % x
            continue
        cur.execute("""SELECT %s_gene.* from %s_gene JOIN %s_%s_orth USING (%s_gene_key)
                                                     JOIN %s_gene USING (%s_gene_key)
                                                     WHERE %s_gene.%s_gene_key = ?""" % (sp1_name,sp2_name,sp1_name,sp2_name,sp2_name,sp1_name,sp1_name,sp2_name,sp2_name),(us_hit[sp2_name+"_gene_key"],))
        x['us_orth'] = cur.fetchone()
        x['us_gene'] = us_hit
        x['us_dist'] = -us_dist

        if not ds_hit:
            print >> sys.stderr, "No downstream gene found for %s" % x
            continue

        cur.execute("""SELECT %s_gene.* from %s_gene JOIN %s_%s_orth USING (%s_gene_key)
                                                     JOIN %s_gene USING (%s_gene_key)
                                                     WHERE %s_gene.%s_gene_key = ?""" % (sp1_name,sp2_name,sp1_name,sp2_name,sp2_name,sp1_name,sp1_name,sp2_name,sp2_name),(ds_hit[sp2_name+"_gene_key"],))
        x['ds_orth'] = cur.fetchone() 
        x['ds_gene'] = ds_hit
        x['ds_dist'] = ds_dist
        try:
            print "us: %s" % (x['us_orth']['name'])
        except:
            pass
        try:
            print "ds: %s" % (x['ds_orth']['name'])
        except:
            pass
                
    overlaps = []

    for y in xls1_d:
        for x in xls2_d:
            if ((('us_gene' in y.keys()) and ('ds_orth' in x.keys())) and
                (y['us_gene'] and x['ds_orth'])):
                if (y['us_gene']['FBid'] == x['ds_orth']['FBid']):
                    overlap_genes = [y['us_gene']]
                    if ((('ds_gene' in y) and y['ds_gene']) and (y['ds_gene'] not in overlap_genes)):
                        overlap_genes.append(y['ds_gene'])
                    if ((('us_orth' in x) and x['us_orth']) and (x['us_orth'] not in overlap_genes)):
                        overlap_genes.append(x['us_orth'])                                           
                    ov_entr = {'overlap_genes' : overlap_genes,
                               'sp1_peaks' : [y],
                               'sp2_peaks' : [x]}
                    k = check_intergenic(overlaps,overlap_genes)
                    if k:
                        if y not in overlaps[k]['sp1_peaks']:
                            overlaps[k]['sp1_peaks'].append(y)
                        if x not in overlaps[k]['sp2_peaks']:
                            overlaps[k]['sp2_peaks'].append(x)
                        for z in overlap_genes:
                            if z['FBid'] not in [x['FBid'] for x in overlaps[k]['overlap_genes']]:
                                print z['FBid']
                                print [x['FBid'] for x in overlaps[k]['overlap_genes']]
                                overlaps[k]['overlap_genes'].append(z)
                    else:
                        overlaps.append(ov_entr)
                    continue
            if ((('ds_gene' in y.keys()) and ('us_orth' in x.keys())) and
                (y['ds_gene'] and x['us_orth'])):
                if (y['ds_gene']['FBid'] == x['us_orth']['FBid']):
                    overlap_genes = [y['ds_gene']]
                    if ((('us_gene' in y) and y['us_gene']) and (y['us_gene'] not in overlap_genes)):
                        overlap_genes.append(y['us_gene'])
                    if ((('ds_orth' in x) and x['ds_orth']) and (x['ds_orth'] not in overlap_genes)):
                        overlap_genes.append(x['ds_orth'])                                           
                    ov_entr = {'overlap_genes' : overlap_genes,
                               'sp1_peaks' : [y],
                               'sp2_peaks' : [x]}
                    k = check_intergenic(overlaps,overlap_genes)
                    if k:
                        if y not in overlaps[k]['sp1_peaks']:
                            overlaps[k]['sp1_peaks'].append(y)
                        if x not in overlaps[k]['sp2_peaks']:
                            overlaps[k]['sp2_peaks'].append(x)
                        for z in overlap_genes: 
                            if z['FBid'] not in [x['FBid'] for x in overlaps[k]['overlap_genes']]:
                                print z['FBid']
                                print [x['FBid'] for x in overlaps[k]['overlap_genes']]
                                overlaps[k]['overlap_genes'].append(z)
                    else:
                        overlaps.append(ov_entr)
                    continue
            if ((('ds_gene' in y.keys()) and ('ds_orth' in x.keys())) and
                (y['ds_gene'] and x['ds_orth'])):
                if (y['ds_gene']['FBid'] == x['ds_orth']['FBid']):
                    overlap_genes = [y['ds_gene']]
                    if ((('us_gene' in y) and y['us_gene']) and (y['us_gene'] not in overlap_genes)):
                        overlap_genes.append(y['us_gene'])
                    if ((('us_orth' in x) and x['us_orth']) and (x['us_orth'] not in overlap_genes)):
                        overlap_genes.append(x['us_orth'])                                           
                    ov_entr = {'overlap_genes' : overlap_genes,
                               'sp1_peaks' : [y],
                               'sp2_peaks' : [x]}
                    k = check_intergenic(overlaps,overlap_genes)
                    if k:
                        if y not in overlaps[k]['sp1_peaks']:
                            overlaps[k]['sp1_peaks'].append(y)
                        if x not in overlaps[k]['sp2_peaks']:
                            overlaps[k]['sp2_peaks'].append(x)
                        for z in overlap_genes:
                            if z['FBid'] not in [x['FBid'] for x in overlaps[k]['overlap_genes']]:
                                print z['FBid']
                                print [x['FBid'] for x in overlaps[k]['overlap_genes']]
                                overlaps[k]['overlap_genes'].append(z)
                    else:
                        overlaps.append(ov_entr)
                    continue
            if ((('us_gene' in y.keys()) and ('us_orth' in x.keys())) and
                (y['us_gene'] and x['us_orth'])):
                if (y['us_gene']['FBid'] == x['us_orth']['FBid']):
                    overlap_genes = [y['us_gene']]
                    if ((('ds_gene' in y) and y['ds_gene']) and (y['ds_gene'] not in overlap_genes)):
                        overlap_genes.append(y['ds_gene'])
                    if ((('ds_orth' in x) and x['ds_orth']) and (x['ds_orth'] not in overlap_genes)):
                        overlap_genes.append(x['ds_orth'])                                           
                    ov_entr = {'overlap_genes' : overlap_genes,
                               'sp1_peaks' : [y],
                               'sp2_peaks' : [x]}
                    k = check_intergenic(overlaps,overlap_genes)
                    if k:
                        if y not in overlaps[k]['sp1_peaks']:
                            overlaps[k]['sp1_peaks'].append(y)
                        if x not in overlaps[k]['sp2_peaks']:
                            overlaps[k]['sp2_peaks'].append(x)
                        for z in overlap_genes:
                            if z['FBid'] not in [x['FBid'] for x in overlaps[k]['overlap_genes']]:
                                print z['FBid']
                                print [x['FBid'] for x in overlaps[k]['overlap_genes']]
                                overlaps[k]['overlap_genes'].append(z)
                    else:
                        overlaps.append(ov_entr)
                    continue
    test = [([(x['name'],x['FBid']) for x in y['overlap_genes']],[(x['chr'],x['loc']) for x in y['sp1_peaks']],[(x['chr'],x['loc']) for x in y['sp2_peaks']]) for y in overlaps]
    
    for x in test:
        print "%s,%s,%s" % x

"""
    for y in xls1_d:
        for x in xls2_d:
            try:
                if (y['us_dist'] < y['ds_dist']) and y['us_gene']['FBid'] == x['ds_orth']['FBid']:
                    print '1'
                    print 'dmel = %s:%s %s => dpse = %s:%s %s %s' % (y['chr'],y['loc'],y['us_dist'],x['chr'],x['loc'],x['ds_dist'],x['ds_orth']['name']) 
                    ov_entr = {'overlap_gene' : (y['us_gene'],x['ds_orth']),
                               'sp1_peaks' : [y],
                               'sp2_peaks' : [x]}
                    if x['us_orth']['FBid'] in overlaps.keys():
                        overlaps[x['ds_orth']['FBid']]['sp1_peaks'].append(y)
                        overlaps[x['ds_orth']['FBid']]['sp2_peaks'].append(x)
                    else:    
                        overlaps[x['ds_orth']['FBid']] = ov_entr
                    continue
            except:
                pass
            try:
                if (y['ds_dist'] < y['us_dist']) and y['ds_gene']['FBid'] == x['us_orth']['FBid']:
                    print '2'
                    print 'dmel = %s:%s %s => dpse = %s:%s %s %s' % (y['chr'],y['loc'],y['ds_dist'],x['chr'],x['loc'],x['us_dist'],x['us_orth']['name']) 
                    ov_entr = {'overlap_gene' : (y['ds_gene'],x['us_orth']),
                               'sp1_peaks' : [y],
                               'sp2_peaks' : [x]}
                    if x['us_orth']['FBid'] in overlaps.keys():
                        overlaps[x['us_orth']['FBid']]['sp1_peaks'].append(y)
                        overlaps[x['us_orth']['FBid']]['sp2_peaks'].append(x)
                    else:    
                        overlaps[x['us_orth']['FBid']] = ov_entr
                    continue
            except:
                pass
            try:
                if (y['ds_dist'] < y['us_dist']) and y['ds_gene']['FBid'] == x['ds_orth']['FBid']:
                    print '3'
                    print 'dmel = %s:%s %s => dpse = %s:%s %s %s' % (y['chr'],y['loc'],y['ds_dist'],x['chr'],x['loc'],x['ds_dist'],x['ds_orth']['name']) 
                    ov_entr = {'overlap_gene' : (y['ds_gene'],x['ds_orth']),
                               'sp1_peaks' : [y],
                               'sp2_peaks' : [x]}
                    if x['us_orth']['FBid'] in overlaps.keys():
                        overlaps[x['ds_orth']['FBid']]['sp1_peaks'].append(y)
                        overlaps[x['ds_orth']['FBid']]['sp2_peaks'].append(x)
                    else:    
                        overlaps[x['ds_orth']['FBid']] = ov_entr
                    continue
            except:
                pass
            try:
                if (y['us_dist'] < y['ds_dist']) and y['us_gene']['FBid'] == x['us_orth']['FBid']:
                    print '4 %s %s' % (y['us_gene']['name'],x['us_orth']['name']) 
                    print 'dmel = %s:%s %s => dpse = %s:%s %s %s' % (y['chr'],y['loc'],y['us_dist'],x['chr'],x['loc'],x['us_dist'],x['us_orth']['name']) 
                    ov_entr = {'overlap_gene' : (y['us_gene'],x['us_orth']),
                               'sp1_peaks' : [y],
                               'sp2_peaks' : [x]}
                    if x['us_orth']['FBid'] in overlaps.keys():
                        print 'a'
                        print ov_entr
                        overlaps[x['us_orth']['FBid']]['sp1_peaks'].append(y)
                        overlaps[x['us_orth']['FBid']]['sp2_peaks'].append(x)
                    else:
                        print 'b'
                        print ov_entr
                        overlaps[x['us_orth']['FBid']] = ov_entr
                    continue
            except:
                pass

    no_overlaps_sp1 = [x for x in xls1_d if (x not in sum(map(lambda y:y['sp1_peaks'],overlaps.values()),[]))]
    no_overlaps_sp2 = [x for x in xls2_d if (x not in sum(map(lambda y:y['sp2_peaks'],overlaps.values()),[]))]

    test = [(x['overlap_gene'][0]['name'],len(x['sp1_peaks']),[(y['chr'],y['loc']) for y in x['sp1_peaks']],len(x['sp2_peaks']),[(y['chr'],y['loc']) for y in x['sp2_peaks']]) for x in overlaps.values()]

    print test
"""


#    print [[(y['chr'],y['loc']) for y in x['sp1_peaks']] for x in overlaps]
#    print [[(y['chr'],y['loc']) for y in x['sp2_peaks']] for x in overlaps]

#    print map(lambda x:(x['chr'],x['loc']),sum(map(lambda y:y['sp2_peaks'],overlaps.values()),[]))
#    print overlaps
#    print no_overlaps_sp1
#    print no_overlaps_sp2


def check_intergenic(overlaps,genes):
    key = None
    for (k,v) in enumerate(overlaps):
        if key:
            break
        for x in [x['FBid'] for x in genes]:
            if x in [y['FBid'] for y in v['overlap_genes']]:
                key = k
    print key
    return key
            

def get_ds_hit(xls_line,gene_hits,cursor):
    min_hit = None
    min_dist = None
    seen = 0
    for x in gene_hits:
        dist = None
        if x['strand'] == "+":
            dist = (x['start'] - xls_line['loc'])
        else:
            dist = (x['end'] - xls_line['loc'])
        if ((not (min_dist)) or (min_dist >= dist)) and dist >= 0:
            min_dist = dist
            min_hit = x
        if x['FBid'] == 'FBgn0073543':
            print xls_line
            print 'dist: %s min_dist: %s' % (dist,min_dist)
        if seen:
            print 'dist: %s min_dist: %s' % (dist,min_dist)

    return (min_hit,min_dist)  

def get_us_hit(xls_line,gene_hits,cursor):
    min_hit = None
    min_dist = None
    seen = 0
    for x in gene_hits:
        dist = None
        if x['strand'] == "+":
            dist = (xls_line['loc'] - x['start'])
        else:
            dist = (xls_line['loc'] - x['end'])
        if ((not (min_dist)) or (min_dist >= dist)) and dist >= 0:
            min_dist = dist
            min_hit = x
        if (x['FBid'] == 'FBgn0073543') or (x['FBid'] == 'FBgn0034430'):
            print xls_line
            print 'dist: %s min_dist: %s %s' % (dist,min_dist,min_hit['name'])
        if seen:
            print 'dist: %s min_dist: %s %s' % (dist,min_dist,min_hit['name'])

    return (min_hit,min_dist)  

if __name__ == "__main__":
    main(sys.argv[1:])
