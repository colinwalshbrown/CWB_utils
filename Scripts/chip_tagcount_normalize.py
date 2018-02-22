#!/Library/Frameworks/Python.framework/Versions/7.1/bin/python
#
# chip_tagcount_normalize.py
# 
# Take a file of mapped reads (BowTie) and create index of tag positions; calculate distribution statistics (n. unique tags, tag frequency distribution); perform read extension and tag normalization
#
#
from __future__ import division

import sys
import re
import gzip as gz
import numpy as np
from optparse import OptionParser
from tables import *
from time import time

# for importing cython modules - TEST
#import pyximport; pyximport.install()
from cython_util_fns import minus_extend,plus_extend

# local modules
from fasta_subseq_2 import *

# globals
FILTERS = Filters(complevel=1)

class BowtieRead(IsDescription):
    def __init__(self,name_size,chr_name_size,readsize):
        self.name = StringCol(name_size)
        self.strand = StringCol(1)
        self.chromosome = StringCol(chr_name_size)
        self.location = UInt32Col()
        self.sequence = StringCol(readsize)
        self.quality_string = StringCol(readsize)

def _main(args):
    # get user options and arguments
    usage = "usage: chip_tagcount_normalize.py [options] <genome_fasta> <bowtie_file>"
    parser = OptionParser(usage=usage)
    parser.add_option("--name","-n")
    parser.add_option("--nfactor","-f",help="normalization factor (int), default 1000000")
    #parser.add_option("--taglimit","-t",help="create separate tag map with no more than limit_n tags per position")
    parser.add_option("--read_format","-r",help="file format for read file - 'sam' or 'bowtie', default = 'sam'")
    parser.add_option("--mapq_min","-q",help="minimum mapq value - use only if format == 'sam'")
    #parser.add_option("--limit_n","-l",help="max number of tags per position in taglimited maps, default = 3")
    parser.add_option("--extend","-e",help="extend tags this number of bases in the 3' direction (~fragment size)")
    #parser.add_option("--tagfrequency","-g",action="store_true",help="output tag frequency distribution information",default=False)
    parser.add_option("--stepsize","-s",help="window size for output sgr file, default 10")
    parser.add_option("--output","-o",action="store_true",help="supress output to wig file(s), only store hdf5",default=True)
    parser.add_option("--store_unique","-u",action="store_true",help="use only unique reads; skips taglimiting",default=False)
    parser.add_option("--zipped_input","-z",action="store_true",help="input file is zipped",default=False)
    parser.add_option("--keep_all_tags","-k",action="store_true",help="keep all tag arrays; otherwise keep only normalized, extended arrays (both lim & non if -l)",default=False)
    options,arguments = parser.parse_args()
    if len(arguments) < 2:
        parser.print_help()
        sys.exit(0)

    # options variables
    fasta = arguments[0]
    bowtie = arguments[1]
    name = (".").join(bowtie.split("/")[-1].split(".")[:-1])
    nfactor = 1000000.0
    taglimited = None
    read_format = "sam"
    #tagfrequency = options.tagfrequency
    store_unique = options.store_unique
    extend = 0
    stepsize = 10
    limit_n = 3
    mapq_min = None
    
    if options.name:
        name = options.name
    if options.nfactor:
        nfactor = float(options.nfactor)
    if options.extend:
        extend = int(options.extend)
    if options.stepsize:
        stepsize = options.stepsize
    #if options.taglimit and (not options.store_unique):
    #    taglimited = True
    #    limit_n = int(options.taglimit)
    if options.read_format:
        read_format = options.read_format
    if (read_format == 'sam') and (options.mapq_min):
        mapq_min = int(options.mapq_min)

    # other variables
    hf5file = openFile(name + ".h5", mode ="w",title=name + " data table")
    chr_arrays = {}
    tags = {}
    limit_tags = {}
    norm_limit_tags = {}
    extend_limit_tags = {}
    norm_tags = {}
    ext_tags = {}
    tdist = []
    nunique = 0
    tag_total = 0
    
    # make arrays to hold tag counts
    chr_arrays,chr_name_size = buildArraysFromFasta(fasta)

    # count tags, calculate tag distribution, and number of unique reads
    tag_arrays,tunique = countTags(bowtie,hf5file,chr_arrays,chr_name_size,read_format,mapq_min,store_unique,options.zipped_input)

    if not options.keep_all_tags:
        hf5file.removeNode(hf5file.root.read_data,recursive=True)

    # sum total mapped tags and normalize to nfactor
    norm_tags = normTags(tag_arrays,nfactor,hf5file,"")

    if not options.keep_all_tags:
        hf5file.removeNode(hf5file.root.tag_counts.raw_tags,recursive=True)

    if extend:
        ext_tags = extendTags(norm_tags,extend,hf5file,"")

    if not options.keep_all_tags:
        hf5file.removeNode(hf5file.root.tag_counts.norm_tags,recursive=True)


"""
    if taglimited:
        print "t"
        limit_tags = limitTags(tag_arrays,limit_n,hf5file)
        norm_limit_tags = normTags(limit_tags,nfactor,hf5file,"Lim")
        if extend:
            extend_limit_tags = extendTags(norm_limit_tags,extend,hf5file,"Lim")

    if not options.keep_all_tags:
        hf5file.removeNode(hf5file.root.read_data,recursive=True)
        hf5file.removeNode(hf5file.root.tag_counts.norm_tags,recursive=True)
        
        hf5file.removeNode(hf5file.root.tag_counts
            
    if tagfrequency:
        tagFrequency(tag_arrays,hf5file,name)

    if tagfrequency:
        

    if output:
        output(name,hf5file,stepsize,tags,"raw_tags")
        output(name,hf5file,stepsize,norm_tags,"norm_tags")
        if extend:
            output(name,hf5file,stepsize,ext_tags,"ext_tags")
        if taglimited:
            output(name,hf5file,stepsize,limit_tags,"limit_tags")
            output(name,hf5file,stepsize,norm_limit_tags,"norm_limit_tags")
            if extend:
                output(name,hf5file,stepsize,extend_limit_tags,"extend_limit_tags")
"""
def buildArraysFromFasta(fasta_name):
    # Create empty arrays, one per chromosome.  Get length information from FastaDB.

    fa = FastaDB()
    fa.openFastaFile(fasta_name)
    chr_name_size = max(map(len,fa.keys()))
    chr_arrays = {}
    for (chr,s) in fa.items():
#        print (chr,s)
        chr_arrays[chr] = {}
        chr_arrays[chr]["+"] = np.zeros(int(s['length']),dtype=np.uint32)
        chr_arrays[chr]["-"] = np.zeros(int(s['length']),dtype=np.uint32)
    return chr_arrays,chr_name_size

#    a_group = hf5file.createGroup(hf5file.root,'raw_counts',"Raw (Not Normalized / Extended) Read Count Arrays")
#    for (chr, in fa.items():
#        hf5file.createArray(a_group,chr,array(np.zeros(),"Chromosome " + chr " Raw Reads")
        
def bowTieReadFactory(name_size,chr_name_size,readsize):
    # Class Factory for making PyTables definition class

    class BowtieRead(IsDescription):
        name = StringCol(name_size)
        strand = StringCol(1)
        chromosome = StringCol(chr_name_size)
        location = UInt32Col()
        sequence = StringCol(readsize)
        quality_string = StringCol(readsize)

    return BowtieRead

def samReadFactory(name_size,chr_name_size,readsize):

        class SamRead(IsDescription):
            name = StringCol(name_size)
            strand = StringCol(1)
            chromosome = StringCol(chr_name_size)
            location = UInt32Col()
            sequence = StringCol(readsize)
            quality_string = StringCol(readsize)
            mapq = UInt32Col()
        return SamRead

def countTags(read_file,hf5file,chr_arrays,chr_name_size,fmt,mapq_min,store_unique,zipped):
    if fmt == "sam":
        return countTagsSam(read_file,hf5file,chr_arrays,chr_name_size,mapq_min,store_unique,zipped)
    elif fmt == "bowtie":
        return countTagsBowtie(read_file,hf5file,chr_arrays,chr_name_size,mapq_min,store_unique,zipped)

def countTagsSam(sam_file,hf5file,chr_arrays,chr_name_size,mapq_min,store_unique,zipped):
    unique_tags = 0
    total_tags = 0
    stored_tags = 0
    mapq_min_tags = 0
    chr_carrays = {}

    sm_file = None
    if zipped:
        sm_file = gz.open(sam_file)
    else:
        sm_file = open(sam_file)
    
    # get field size information from first line of BowTie file
    first_line = sm_file.readline().split("\t")
#    print first_line
    name_size = len(first_line[0])
    readsize = len(first_line[9])
    sam_class = samReadFactory(name_size,chr_name_size,readsize)
    
    # create pytable for BowTie info
    group = hf5file.createGroup("/",'read_data','Raw Reads')
    table = hf5file.createTable(group,'sam',sam_class,'SAM Reads',filters=FILTERS)
    sm_file.close()

    r = table.row

    sm_file = None
    if zipped:
        sm_file = gz.open(sam_file)
    else:
        sm_file = open(sam_file)

    for l in sm_file:
        if re.match("^@",l):
            continue
        line = l[:-1].split("\t")
#        print line
        chrom = line[2]
        location = int(line[3])
        strand = None
        if int(line[1]) | 16 == int(line[1]):
            strand = "-"
        else:
            strand = "+"

        r['name'] = line[0]
        r['strand'] = strand
        r['chromosome'] = chrom
        r['location'] = location
        r['sequence'] = line[9]
        r['quality_string'] = line[10]
        r['mapq'] = int(line[4])

        total_tags += 1
        if total_tags % 100000 == 0:
            print (total_tags,unique_tags)
        if ((chr_arrays[chrom][strand][location]) == 0):
            unique_tags += 1
            if mapq_min and (r['mapq'] >= mapq_min):
                r.append()
                chr_arrays[chrom][strand][location] += 1
                stored_tags += 1
            elif mapq_min:
                mapq_min_tags += 1
            elif not mapq_min:
                r.append()
                chr_arrays[chrom][strand][location] += 1
                stored_tags += 1
        elif not (store_unique):
            if mapq_min and (r['mapq'] >= mapq_min):
                r.append()
                chr_arrays[chrom][strand][location] += 1
                stored_tags += 1
            elif mapq_min:
                mapq_min_tags += 1
            elif not mapq_min:
                r.append()
                chr_arrays[chrom][strand][location] += 1
                stored_tags += 1
        
                
#        print chr_arrays[line[2]][int(line[3])-50:int(line[3])+50]

    print "*** final tag totals: %d tags, %d unique, %d suppressed by mapq_min, %d stored" % (total_tags,unique_tags,mapq_min_tags,stored_tags)
    table.flush()

    # store arrays in pytable
    tag_group = hf5file.createGroup("/",'tag_counts','Tag Counts (by chromosome)')
    raw_group = hf5file.createGroup(tag_group,'raw_tags','Raw (non-normalized, non-extended) Tag Count')

    for chr,a in chr_arrays.items():
        print "storing %s" % (chr,)

        chr_carrays[chr] = {"-":None,"+":None}

        chr_carrays[chr]["+"] = hf5file.createCArray(raw_group,"chr"+chr+"_plus",UInt32Atom(),(len(a['+']),),title=chr + "Plus Strand",filters=FILTERS)
        chr_carrays[chr]["-"] = hf5file.createCArray(raw_group,"chr"+chr+"_minus",UInt32Atom(),(len(a['-']),),title=chr + "Minus Strand",filters=FILTERS)

        chr_carrays[chr]["+"][0:] = a['+']
        chr_carrays[chr]["-"][0:] = a['-']

    return chr_carrays,unique_tags

        
def countTagsBowtie(bowtie_file,hf5file,chr_arrays,chr_name_size,mapq_min,store_unique,zipped):
    unique_tags = 0
    total_tags = 0
    stored_tags = 0
    chr_carrays = {}

    bt_file = None
    if zipped:
        bt_file = gz.open(bowtie_file)
    else:
        bt_file = open(bowtie_file)
    
    # get field size information from first line of BowTie file
    first_line = bt_file.readline().split("\t")
    print first_line
#    print first_line
    name_size = len(first_line[0])
    readsize = len(first_line[4])
    bt_class = bowTieReadFactory(name_size,chr_name_size,readsize)
    
    # create pytable for BowTie info
    group = hf5file.createGroup("/",'read_data','Raw Reads')
    table = hf5file.createTable(group,'bowtie',bt_class,'Bowtie Reads',filters=FILTERS)
    bt_file.close()

    r = table.row

    bt_file = None
    if zipped:
        bt_file = gz.open(bowtie_file)
    else:
        bt_file = open(bowtie_file)

    for l in bt_file:
        line = l[:-1].split("\t")
#        print line

        location = int(line[3])
        chrom = line[2]
        strand = line[1]

        r['name'] = line[0]
        r['strand'] = line[1]
        r['chromosome'] = line[2]
        r['location'] = location
        r['sequence'] = line[4]
        r['quality_string'] = line[5]
        #r.append()

        if ((chr_arrays[chrom][strand][location]) == 0):
            unique_tags += 1
            r.append()
            chr_arrays[chrom][strand][location] += 1
            stored_tags += 1
        elif not (store_unique):
            r.append()
            chr_arrays[chrom][strand][location] += 1
            stored_tags += 1
            
        total_tags += 1
        if total_tags % 100000 == 0:
            print (total_tags,unique_tags)
#        print chr_arrays[line[2]][int(line[3])-50:int(line[3])+50]

    print "*** final tag totals: %d tags, %d unique, %d stored" % (total_tags,unique_tags,stored_tags)
    table.flush()

    # store arrays in pytable
    tag_group = hf5file.createGroup("/",'tag_counts','Tag Counts (by chromosome)')
    raw_group = hf5file.createGroup(tag_group,'raw_tags','Raw (non-normalized, non-extended) Tag Count')

    for chr,a in chr_arrays.items():
        print "storing %s" % (chr,)

        chr_carrays[chr] = {"-":None,"+":None}

        chr_carrays[chr]["+"] = hf5file.createCArray(raw_group,"chr"+chr+"_plus",UInt32Atom(),(len(a['+']),),title=chr + "Plus Strand",filters=FILTERS)
        chr_carrays[chr]["-"] = hf5file.createCArray(raw_group,"chr"+chr+"_minus",UInt32Atom(),(len(a['-']),),title=chr + "Minus Strand",filters=FILTERS)

        chr_carrays[chr]["+"][0:] = a['+']
        chr_carrays[chr]["-"][0:] = a['-']

    return chr_carrays,unique_tags

def normTags(tags,nfactor,hf5file,label):

    norm_group = hf5file.createGroup(hf5file.root.tag_counts,'norm%s_tags' % (label,),'Normalized %s Tags' % (label,))

    norm_tags = {}
    tag_total = 0
    norm_tag_count = 0

    # Calculate Tag Total
    for (chr,a_obj) in tags.items():     
        
        print "counting %s" % (chr,)
        
        plus = a_obj["+"].read()
        minus = a_obj["-"].read()

        tag_total += plus.sum()
        tag_total += minus.sum()
        
        #print "unread + sum: %f" % (a_obj["+"].sum())
        print "total: %f" % (tag_total,)

    # Apply normalization chr by chr and store
    for (chr,a_obj) in tags.items():
        plus = a_obj["+"].read()
        minus = a_obj["-"].read()
        norm_tags[chr] = {}
        norm_minus = minus * (nfactor / tag_total)
        norm_plus = plus * (nfactor / tag_total)
#        norm_tags[chr]['-'] = hf5file.createArray(norm_group,"chr"+chr + "_minus", np.array(norm_minus), chr + " Minus strand") 
#        norm_tags[chr]['+'] = hf5file.createArray(norm_group,"chr"+chr + "_plus", np.array(norm_plus), chr + " Plus strand")
        norm_tags[chr]['-'] = hf5file.createCArray(norm_group,"chr"+chr + "_minus",Float64Atom(), (len(norm_minus),), title=chr + " Minus strand",filters=FILTERS) 
        norm_tags[chr]['+'] = hf5file.createCArray(norm_group,"chr"+chr + "_plus",Float64Atom(), (len(norm_plus),), title=chr + " Plus strand",filters=FILTERS)
        norm_tags[chr]['-'][0:] = norm_minus
        norm_tags[chr]['+'][0:] = norm_plus
        print "norm tags for %s: %f" % (chr,(norm_minus.sum() + norm_plus.sum()))
        norm_tag_count += (norm_minus.sum() + norm_plus.sum())

    print "Starting Tag Total: %d" % (tag_total,)
    print "Normalized Tag Total: %f" % (norm_tag_count,)
    return norm_tags

def extendTags(tags,ext,hf5file,label):
    
    ext_group = hf5file.createGroup(hf5file.root.tag_counts,'ext%s_tags' % (label,),'Extended, Normalized %s Tags' % (label,))
    ext_tags = {}

    for (chr,a_obj) in tags.items():
        print "extending %s" % (chr,)

        # copy normalized reads into memory
        plus = a_obj["+"].read()
        minus = a_obj["-"].read()

#        ext_a_c = np.zeros(len(plus))
#        ext_a_n = np.zeros(len(plus))
        a_len = len(plus)

 #       print "start cython"
        c_st = time()
        # Extend Strands (external cython methods)
        print "1"
        ext_a_c_m = minus_extend(minus,a_len,ext)
 #       print "c minus sum: %f" % ext_a_c_m.sum()
        print "2"
        ext_a_c_p = plus_extend(plus,a_len,ext)
 #       print "c plus sum: %f" % ext_a_c_p.sum()
        print "3"
        ext_a_c_m += ext_a_c_p
        print "done extending - %f elapsed" % (time() - c_st)

        # Note that this is no longer split by strand
        ext_tags[chr] = hf5file.createCArray(ext_group,"ext%s_%s" % (label,chr),Float64Atom(),(len(ext_a_c_m),),title="Extended, Normalized %s Tags" % (label,),filters=FILTERS)
        ext_tags[chr][0:] = ext_a_c_m

    return ext_tags

def limitTags(tags,limit,hf5file):
    
    lim_group = hf5file.createGroup(hf5file.root.tag_counts,"limit_tags","Limited Tags (to 'limit' tags / position)")
    lim_tags = {}

    for (chr,a_obj) in tags.items():

        # copy raw read arrays into memory
        plus = a_obj["+"].read()
        minus = a_obj["-"].read()

        # clip arrays 
        print "limiting %s - %s tags, %d tag limit" % (chr,(plus.sum() + minus.sum()),limit)
        lim_a_p = plus.clip(min=0.0,max=limit)
        lim_a_m = minus.clip(min=0.0,max=limit)
        print "Done - %s Tags after limiting" % (lim_a_p.sum() + lim_a_m.sum())

        # Write chromosome array back to hf5file
        lim_tags[chr] = {}
#        lim_tags[chr]["+"] = hf5file.createArray(lim_group,"lim_%s_plus" % (chr,),np.array(lim_a_p),"Limited Raw %s Tags; Plus Strand" % (chr,))
#        lim_tags[chr]["-"] = hf5file.createArray(lim_group,"lim_%s_minus" % (chr,),np.array(lim_a_m),"Limited Raw %s Tags; Minus Strand" % (chr,))
        lim_tags[chr]["+"] = hf5file.createCArray(lim_group,"lim_%s_plus" % (chr,),UInt32Atom(),(len(lim_a_p),),title="Limited Raw %s Tags; Plus Strand" % (chr,),filters=FILTERS)
        lim_tags[chr]["-"] = hf5file.createCArray(lim_group,"lim_%s_minus" % (chr,),UInt32Atom(),(len(lim_a_m),),title="Limited Raw %s Tags; Minus Strand" % (chr,),filters=FILTERS)
        lim_tags[chr]["+"][0:] = lim_a_p
        lim_tags[chr]["-"][0:] = lim_a_m
    
    return lim_tags

def tagFrequency(tags,hf5file,name):

    hist_group = hf5file.createGroup(hf5file.root,"histograms","Histograms of tag frequency counts (for non-normalized data)")
    
    h_file = name + "_hist.txt"
    hists = {}
    h_max = 0

    # find genome-wide tag max
    print "getting dataset max...."
    for (chr,a_obj) in tags.items():
        
        print chr
        plus = a_obj["+"].read()
        minus = a_obj["-"].read()
        
        ch_max = max(plus.max(),minus.max())
        print "chr %s: %d h_max: %d" % (chr,ch_max,h_max)
        if ch_max > h_max:
            print "C: %f H: %f chr: %s" % (ch_max,h_max,chr)
        h_max = int(max(h_max,ch_max))
    print "done - max == %d"  % (h_max,)
        
    h_tot = np.zeros(h_max,dtype=np.uint64)
    bins = np.arange(h_max + 1)

    for (chr,a_obj) in tags.items():
        
        print "Calculating Histogram for %s" % (chr,)
        # copy raw read arrays into memory
        plus = a_obj["+"].read()
        minus = a_obj["-"].read()

        # get max lengths for plus and minus
        #p_max = plus.max()
        #m_max = minus.max()
        
        #h_mx = max(len(tot_hist),p_max,m_max)
        #bins = np.arange(h_max)

        pl_h = np.histogram(plus,bins)[0]
        min_h = np.histogram(minus,bins)[0]

        combo_hist = pl_h + min_h
        
        # zeros are counting of zeros on p and m strands
        h_tag_tot = pl_h.sum()
        h_zeros = pl_h.sum() - (sum(pl_h[1:]) + sum(min_h[1:]))
        combo_hist[0] = h_zeros
        
        hists[chr] = combo_hist
        print "c:",
        print combo_hist
        print "h_tot:",
        print h_tot
        h_tot += combo_hist
        print "Done"
    
    h_out = open(h_file,"w")
    h_count = h_tot.sum()
    print "Hist total tags: %d" % (h_count,)
    h_norm = h_tot / h_count
    print "Norm hist total: %f" % (h_norm.sum())
    raw_hist = hf5file.createArray(hist_group,"raw_counts",np.array(h_tot),"Histogram of raw tag counts for non-normalized data")
    norm_hist = hf5file.createArray(hist_group,"norm_counts",np.array(h_norm),"Histogram of raw tag counts, normalized by total so sum == 1")
    
    print >> h_out, "\t".join(map(str,range(0,len(h_tot))))
    print >> h_out, "\t".join(map(str,h_tot))
    print >> h_out, "\t".join(map(str,h_norm))
    
    return (raw_hist,norm_hist)
        
        
        
                    
        

#        print "start native"
#        n_st = time()
#       ext_a_n = minus_ext_native(minus,ext_a_n,a_len,ext)
#        ext_a_n = plus_ext_native(plus,ext_a_n,a_len,ext)
#        print "done - %f elapsed" % (time() - n_st)

#        diff_a = ext_a_c - ext_a_n
#        print "Cython sum: %d Native Sum: %d" % (ext_a_c.sum(),ext_a_n.sum())
#        print "Diff Sum: %f" % (diff_a.sum())

#            if (x % 100000 == 0) and (minus[x] > 0):
#                print minus[lower:x]
#                print ext_a[lower:x]
        
#            if (x % 100000 == 0):
#                print plus[(x+1):upper]
#                print ext_a[(x+1):upper]

# native extension methods - way slow!
"""
def minus_ext_native(minus,ext_a,a_len,ext):
    for x in range(a_len):
        lower = x - ext
        if lower < 0:
            lower = 0
        ext_a[lower:x] += minus[x]
    return ext_a

def plus_ext_native(plus,ext_a,a_len,ext):
    for x in range(a_len):
        upper = x + ext + 1
        if upper > a_len:
            upper = a_len
        ext_a[(x+1):upper] += plus[x]
    return ext_a
"""

if __name__ == "__main__":
    _main(sys.argv[1:])
