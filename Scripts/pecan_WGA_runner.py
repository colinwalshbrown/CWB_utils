#!/usr/bin/env python

import sys
import alignment
import fasta_subseq_2
import multiprocessing as mp
import subprocess as sub
import itertools as it
import copy as cp
import re
import os
from Bio import Phylo
from cStringIO import StringIO

PECANPATH = "/Users/cwbrown/src/pecan_v0.8/pecan_v0.8.jar"

def _main(args):
    
    usage = "pecan_WGA_runner.py <genome_fastas_file> <treefile> <mercator_map> <mercator_genome_order> <outdir>"
    if len(args) != 5:
        print usage
        sys.exit(0)

    genome_dict = {}
    genome_order = []
    mercator_genome_order = []
    map_lines = []

    mercator_order_file = open(args[3])
    mercator_genome_order = mercator_order_file.readline().split()

    for ln in open(args[0]):
        (sp,fasta) = ln[:-1].split()
        genome_order.append(sp)
        genome_dict[sp] = fasta_subseq_2.FastaDB()
        genome_dict[sp].openFastaFile(fasta)
    
    tree_dict = generate_trees(genome_dict,args[1])

    for ln in open(args[2]):
        map_ln = ln[:-1].split()
        (species,map_dict) = make_map_dict(map_ln[1:],genome_order,mercator_genome_order)
        tree_obj = tree_dict[species]
        tree_strIO = StringIO()
        Phylo.write(tree_obj,tree_strIO,"newick")
        fastas = [(sp,genome_dict[sp]) for sp in genome_order if (sp in species)]
        map_entry = {'map_dict' : map_dict,
                     'tree' : tree_strIO.getvalue(),
                     'map_idx' : int(map_ln[0]),
                     'fastas' : fastas}
        map_lines.append(map_entry)

    os.chdir(args[4])

    pool = mp.Pool(2)

    #pool.map(run_aln_mapline, map_lines)
    pool.map(run_aln_mapline, map_lines)

    print "ALL DONE!"

def generate_trees(genome_dict, treefile):

    spp = genome_dict.keys()

    tree_comb = [tuple(sorted(spp))]
    tree_dict = {}

    #tree_file = open(treefile)
    #tree_string = tree_file.readline()
    #tree_string = "((ms11, (pid18, dgi2)), (SK-93-1035, ((f62, fa1090), (SK-92-679, (1291, (fa6140, (pid24-1, dgi18)))))));\n"
    #print tree_string[:-2]
    base_tree = Phylo.read(treefile,"newick")
    #print "Base: "
    #print base_tree
    #print "================="
    
    for i in range(1,len(sorted(genome_dict.keys()))):
        tree_comb.extend(it.combinations(sorted(genome_dict.keys()),i))

    for com in tree_comb:
        del_sp = [x for x in spp if x not in com]
        #print del_sp
        pruned = cp.deepcopy(base_tree)
        for s in del_sp:
             pruned.prune(s)
        tree_dict[com] = pruned   

    return tree_dict
    

def make_map_dict(map_line,genome_order,mercator_genome_order):

    species = []
    sp_loc_dict = {}
    ix = 0
    #print map_line
    for sp in mercator_genome_order:

        next = ix + 4
        loc_str = []
        if (next) >= len(map_line):
            loc_str = map_line[ix:]
        else:
            loc_str = map_line[ix:ix+4]
        
        if loc_str[0] == "NA":
            ix = next
            continue

        else:
            loc_dict = {'chr':loc_str[0],
                        'start':int(loc_str[1]),
                        'end':int(loc_str[2]),
                        'strand':loc_str[3]}
            #if loc_dict['start'] <= 0:
            #    loc_dict['start'] == 1
            sp_loc_dict[sp] = loc_dict
            species.append(sp)
            ix = next

    return (tuple(sorted(species)),sp_loc_dict)
            

def run_aln_mapline(map_line):

    files = []
    for (sp,fasta) in map_line['fastas']:
        #print >> sys.stderr, (map_line['map_idx'],'a')
        sp_map = map_line['map_dict'][sp]
        #print >> sys.stderr, (map_line['map_idx'],'b')
        outname = str(map_line['map_idx']) + "_" + sp + ".fa"
        out = open(outname,"w")
        #print >> sys.stderr, (map_line['map_idx'],'c')
        print >> sys.stderr, sp_map
        print >> sys.stderr, sp
        seq = fasta[sp_map['chr']][sp_map['start']:sp_map['end'] - 1]
        #print >> sys.stderr, (map_line['map_idx'],'d')
        if sp_map['strand'] == "-":
            print sp + " revcomp"
            seq = fasta_subseq_2.revcomp(seq)
        #print >> sys.stderr, (map_line['map_idx'],'e')
        print >> out, ">%s %s-%d:%d" % (sp,sp_map['chr'],sp_map['start'],sp_map['end'])
        print >> out, seq
        print "Wrote fasta %s" % (outname,)
        #print (map_line['map_idx'],'f')
        out.close()
        files.append(outname)
        #print >> sys.stderr, (map_line['map_idx'],'g')
        

    call_pecan = ["java","-classpath",PECANPATH,"-Xmx2000m","bp.pecan.Pecan","-E",map_line['tree'] + ";","-F"]
    call_pecan.extend(files)
    call_pecan.extend(["-G",str(map_line['map_idx']) + ".mfa"])
    print "Running pecan with command: %s" % (" ".join(call_pecan))
    print "Starting alignment %d..." % (map_line['map_idx'],)
    sub.check_call(call_pecan)
    print "Alignment %d finished" % (map_line['map_idx'],)
    map(os.remove,files)

if __name__ == "__main__":
    _main(sys.argv[1:])
