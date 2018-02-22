#!/usr/bin/env python

import sys
import fasta_subseq
import re
import sqlite3
import matplotlib.pyplot as plt

def main(args):
    
    if len(args) <4 :
        print "motif_coocurrance.py: <xls> <patser_db> <order_str=0,1,0...> <window> <genome> <outbase> <motif1> <motif2> [...]"
        sys.exit(1)

    conn_sp1 = sqlite3.connect(args[1])
    conn_sp1.row_factory = sqlite3.Row
    cur_sp1 = conn_sp1.cursor()
    peaks = args[0]
    window = int(args[3])
    outbase = args[5]
    mtx_names = args[6:]
    print mtx_names

    genome = fasta_subseq.FastaDB()
    genome.openFastaFile(args[4])

    order = args[2].split(",")

    matrices = []
    pks = []

    for mtx in mtx_names:
        cur_sp1.execute("""SELECT matrix_key from matrix where name = ?""",(mtx,))
        matrices.append(cur_sp1.fetchall()[0][0])

    for ln in open(peaks):
        if re.search("#",ln):
            continue
        ln_spl = ln[:-1].split("\t")
        for (idx,pk) in enumerate(ln_spl[4][:-1].split(",")):
            if not pk:
                continue
            pk_loc ={'chr':ln_spl[1],'reg_id':ln_spl[0],'loc':int(pk),'hits':{},'start':int(pk) - window/2,'end':int(pk)+window/2}
            print({'chr':ln_spl[1],'reg_id':ln_spl[0],'loc':int(pk)})
            for mtx in matrices:
                cur_sp1.execute("""SELECT * from patser_hit where ((chr = ?) and (start between ? and ?) and (matrix_key = ?))""",(pk_loc['chr'],pk_loc['start'],pk_loc['end'],mtx))
                hits = cur_sp1.fetchall()
                if len(hits) > 0:
                    pk_loc['hits'][mtx] = hits
    
            pks.append(pk_loc)
            print pk_loc
    
    #pks.sort(key=(lambda x:len(x['hits'])))

    hits_5_4_3 = []
    hits_2 = []
    hits_1 = []

    for (idx,x) in enumerate(pks):
        ref_mtx = None
        orient = None
        hits = None
        print len(x['hits'])
        if len(x['hits']) > 2:
            hits = hits_5_4_3
        elif len(x['hits']) == 2:
            hits = hits_2
        elif len(x['hits']) == 1:
            hits = hits_1
        else:
            continue
        for y in range(0,len(order)):
            print matrices[y]
            print x['hits'].keys()
            if matrices[y] in x['hits'].keys():
                ref_mtx = matrices[y]
                orient = order[y]
                break
        if not ref_mtx or not orient:
            print "no hits for %s" % (str(x),)
            continue
        #hits['peaks'] += 1
        print "ref_mtx: %s orient: %s" % (ref_mtx,orient)
        sort_pks = sorted(x['hits'][ref_mtx],key=(lambda z:abs(z['start'] - x['loc'])))
        print sort_pks
        middle_pk = sort_pks[0]
        print "start: %s:%s end %s:%s" % (x['chr'],x['start'],x['chr'],x['end'])
        peak = []
        if (middle_pk['strand'] == orient):
            x['strand'] = '+'
            for (name,hts) in x['hits'].items():
                for h in hts:
                    end_dist = h['end'] - x['start']
                    start_dist = h['start'] - x['start']
                    hit_mid = (start_dist + end_dist) / 2
                    strand = h['strand']
#                    print (h['start'],h['end'])
                    flp_h = {'mtx':mtx_names[matrices.index(name)],'chr':h['chr'],'start':h['start'],'end':h['end'],'mid':hit_mid,'start_dist':start_dist,'end_dist':end_dist,'strand':h['strand'],'orig_strand':h['strand']}
#                    print flp_h
#                    hits['peaks'] += 1
#                    print hits['peaks']
#                    if (mtx_names[matrices.index(name)],strand) in hits.keys():
#                        hits[(mtx_names[matrices.index(name)],strand)].append((hits['peaks'],hit_mid))
#                    else:
#                        hits[(mtx_names[matrices.index(name)],strand)] = [(hits['peaks'],hit_mid)]
#                    if (mtx_names[matrices.index(name)],strand) in hits.keys():
#                        hits[(mtx_names[matrices.index(name)],strand)].append((hits['peaks'],hit_mid))
#                    else:
#                        hits[(mtx_names[matrices.index(name)],strand)] = [(hits['peaks'],hit_mid)]
                    if h == middle_pk:
                        peak.insert(0,flp_h)
                    else:
                        peak.append(flp_h)
            hits.append({'hits':peak,'id':x['reg_id'],'chr':x['chr'],'start':x['start'],'end':x['end'],'strand':x['strand']})
        else:
            flipped_hits = []
            x['strand'] = '-'
            for (name,hts) in x['hits'].items():
                for h in hts:
                    print "flip:"
                    print h
                    end_dist = (x['end'] - x['start']) - (h['start'] - x['start'])
                    st_dist = x['end'] - h['end']
                    hit_mid = (st_dist + end_dist)/2
                    strand = None
                    if (h['strand'] == "+"):
                        strand = "-"
                    else:
                        strand = "+"
                    flp_h = {'mtx':mtx_names[matrices.index(name)],'chr':h['chr'],'mid':hit_mid,'start':h['start'],'end':h['end'],'start_dist':st_dist,'end_dist':end_dist,'orig_strand':h['strand'],'strand':strand}
 #                   print flp_h
 #                    hits['peaks'] += 1
 #                   print hits['peaks']
                    if h == middle_pk:
                        peak.insert(0,flp_h)
                    else:
                        peak.append(flp_h)
            hits.append({'hits':peak,'id':x['reg_id'],'chr':x['chr'],'start':x['start'],'end':x['end'],'strand':x['strand']})
            
        print "hits: %s" % (hits,)

    

    plot_chars = {(mtx_names[0],"+"):"m>",
                  (mtx_names[0],"-"):"m<",
                  (mtx_names[1],"+"):"r>",
                  (mtx_names[1],"-"):"r<",
                  (mtx_names[2],"+"):"g>",
                  (mtx_names[2],"-"):"g<",
                  (mtx_names[3],"+"):"b>",
                  (mtx_names[3],"-"):"b<",
                  (mtx_names[4],"+"):"y>",
                  (mtx_names[4],"-"):"y<"}

    hits_5_4_3.sort(key=(lambda x:x['hits'][0]['mid']))
    hits_2.sort(key=(lambda x:x['hits'][0]['mid']))
    hits_1.sort(key=(lambda x:x['hits'][0]['mid']))

    print hits_5_4_3
    print hits_2
    print hits_1

    motif_1_gt3 = open(outbase + "_motif_1_gt3-others.txt","w")
    motif_1_2 = open(outbase + "_motif_1_2-others.txt","w")
    motif_1_1 = open(outbase + "_motif_1_no-others.txt","w")

    seq_1_gt3 = open(outbase + "_motif_1_gt3-others_seqs.txt","w")
    seq_1_2 = open(outbase + "_motif_1_2-others_seqs.txt","w")
    seq_1_1 = open(outbase + "_motif_1_no-others_seqs.txt","w")

    for (ht,mot,seqout) in ((hits_5_4_3,motif_1_gt3,seq_1_gt3),(hits_2,motif_1_2,seq_1_2),(hits_1,motif_1_1,seq_1_1)):
        for pk in ht:
            reg_seq = genome[pk['chr']]['sequence'][pk['start']:pk['end']]
            if pk['strand'] == "-":
                reg_seq = fasta_subseq.revcomp(reg_seq)
            print >> seqout, ">%s:%s-%s:%s" % (pk['chr'],pk['start'],pk['end'],pk['strand'])
            print >> seqout, reg_seq
            for h in pk['hits']:
                if h['mtx'] == mtx_names[0]:
                    seq = genome[h['chr']]['sequence'][h['start']:h['end']]
                    if h['orig_strand'] == "-":
                        seq = fasta_subseq.revcomp(seq)
                    print >> mot,seq
                    

    plt.figure(1,figsize=(10,14))

    plot_set = ([x['hits'] for x in hits_5_4_3])

    lists = {}
    labels = [x['id'] for x in hits_5_4_3] 

    for (idx,hits) in enumerate(plot_set):
        if len(hits) == 0:
            continue
        
        for x in mtx_names:
            plus_hts = []
            minus_hts = []
            for (i,pk) in enumerate(hits):
                print pk
                plus_hts.extend([(idx,z['mid']) for z in (pk,) if (z['mtx'] == x and z['strand'] == "+")])
                minus_hts.extend([(idx,z['mid']) for z in (pk,) if (z['mtx'] == x and z['strand'] == "-")])
            if (x,"+") in lists.keys():
                lists[(x,"+")].extend(plus_hts)
            else:
                lists[(x,"+")] = plus_hts
            if (x,"-") in lists.keys():
                lists[(x,"-")].extend(minus_hts)
            else:
                lists[(x,"-")] = minus_hts
            print plus_hts
            print minus_hts
            print lists.keys()
            print lists.items()

    for (k,dr) in lists.items():
        print "hits: %s %s" % (k,dr)
        if len(dr) < 1:
            continue
        y,x=zip(*dr)
        print k
        print x
        print y
#                plt.subplot('1%s%s' % (len(plot_set),idx + 1))
        plt.plot(x,y,plot_chars[k],markersize=10.0,label=k[0])
    plt.yticks(range(len(hits_5_4_3)),labels)
    plt.xlim((-20,window+20))
    plt.grid(True,markevery=1)

    plt.savefig(outbase + "5_4_3_hits.svg",format="svg")

    print "3-5_hits: %s 2_hits: %s 1_hits: %s" % (len(hits_5_4_3),len(hits_2),len(hits_1))
    plt.show()

    

if __name__ == "__main__":
    main(sys.argv[1:])
    
