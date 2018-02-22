#!/usr/bin/env python
######################
#
#   chipseq_processing.py:
#   Module w/ classes and methods for chipseq processing pipeline
#   Uses pytables to keep track of read counts
#   Numpy / cython for doing base-by-base operations
#

from __future__ import division

import sys
import re
import numpy as np
import datetime
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
SPECIES_FASTAS = {'dpse' : "~/Data/genomes/dpse-current.fa",
                  'dmel' : "~/Data/genomes/dmel-current.fa"}
TEMP_NAME = "TEMP_SeqAnalysisTables_%s.hf5" % (datetime.date.today(),)

######################
#
#   Classes for sequencing data 
#

class SeqAnalysisTable():
    def __init__(self,hf5=None):
        # create table in memory by default
        self.hf5_file = hf5
        self.delete_on_close = True
        self.hf5_object = None
        self.array_group = None
        self.read_group = None
        self.read_tables = {}
        self.count_arrays = {}

        if not self.hf5_file:
            self.hf5_file = TEMP_NAME
            self.hf5_object = self.open_hf5()
            # create groups for reads and count arrays
            self.read_group = self.hf5_object.createGroup("/",'read_group','Raw Reads')
            self.array_group = self.hf5_object.createGroup("/",'array_group','Tag Count Arrays')
        else:
            self.delete_on_close = False
            self.hf5_object = self.open_hf5(self)
            self.read_group = self.hf5_object.read_group
            self.array_group = self.hf5_object.array_group
            self.count_arrays = self.get_existing_arrays()
            self.read_tables = self.get_existing_read_tables()

    def get_existing_arrays(self):
        pass
        
    def open_hf5(self):
        return openFile(self.hf5_file, mode ="w",title="Sequencing data table")

    def open_reads(self,readfile,readformat,chr_name_size=32):
        readobj = None
        
        if readformat.lower() == "sam":
            readobj = self.open_sam_reads(readfile,chr_name_size)
        elif readformat.lower() == "bowtie":
            readobj = self.open_bowtie_reads(readfile,chr_name_size)

        return readobj

    def open_sam_reads(self,readfile,chr_name_size):

        reads = open(readfile)
    
        # get field size information from first line of BowTie file
        first_line = reads.readline().split("\t")
        #    print first_line
        name_size = len(first_line[0])
        readsize = len(first_line[9])
        read_class = self.samReadFactory(name_size,chr_name_size,readsize)

        read_coll = ReadCollection(self,read_class)
        reads.close()

        for l in open(readfile):
            if not re.match("^@",l):
                read_coll.add_read(l)

        print "1"
        read_coll.read_table.flush()

        idx = str(read_coll.read_tbl_idx)
        self.read_tables["reads_" + idx] = read_coll
        print read_coll.items()
        
        return read_coll

    def build_tag_array(self,read_table,species=None,genome_fasta=None,tag_limit=None,name=None):

        tag_array = GenomeTrack(self,species=species,genome_fasta=genome_fasta,tag_limit=tag_limit,name=name)

        for r in read_table:
            tag_array.add_tag(r['chromosome'],r['location'],r['strand'])

    def __del__(self):
        self.hf5_object.close()
        if self.delete_on_close:
            os.remove(self.hf5_file)            

    def list_reads(self):
        return self.read_tables.values()

    def list_arrays(self):
        return self.read_tables.values()

    def get_reads(self,key):
        return self.read_tables[key]

    def get_array(self,key):
        return self.array_tables[key]

    def samReadFactory(self,name_size,chr_name_size,readsize):

        class SamRead(IsDescription):
            name = StringCol(name_size)
            strand = StringCol(1)
            chromosome = StringCol(chr_name_size)
            location = UInt32Col()
            sequence = StringCol(readsize)
            quality_string = StringCol(readsize)
            mapq = UInt32Col()

        return SamRead

    def bowTieReadFactory(self,name_size,chr_name_size,readsize):
        # Class Factory for making PyTables definition class for Bowtie lines

        class BowtieRead(IsDescription):
            name = StringCol(name_size)
            strand = StringCol(1)
            chromosome = StringCol(chr_name_size)
            location = UInt32Col()
            sequence = StringCol(readsize)
            quality_string = StringCol(readsize)

        return BowtieRead

class ReadCollection():
    def __init__(self,parent_table,read_class,name=None,table=None):
        self.parent_table = parent_table
        self.read_tbl_idx = len(parent_table.list_reads())
        self.read_class = read_class
        self.read_class_name = read_class.__name__
        self.read_total = 0
        self.name = None
        if not name:
            self.name = "reads_" + str(self.read_tbl_idx)
        else:
            self.name = name
        self.read_table = None
        if not table:
            self.read_table = self.parent_table.hf5_object.createTable(self.parent_table.read_group,self.name,read_class,'Raw Reads ' + self.name)
        else:
            self.read_table = table
 
    def add_read(self,string):
        if self.read_class_name == 'BowtieRead':
            self.parse_bowtie_read(string)
        elif self.read_class_name == 'SamRead':
            self.parse_sam_read(string)
        self.read_total += 1
        if (self.read_total % 1000 == 0):
            print "processed %d reads, stored %d" % (self.read_total,len(self.read_table))
            print string
        if (self.read_total < len(self.read_table)):
            self.read_total = len(self.read_table)
            print string

    def parse_sam_read(self,string):
        line = string[:-1].split("\t")
        r = self.read_table.row

        strand = None
        if int(line[1]) | 16 == int(line[1]):
            strand = "-"
        else:
            strand = "+"

        r['name'] = line[0]
        r['strand'] = strand
        r['chromosome'] = line[2]
        r['location'] = int(line[3])
        r['sequence'] = line[9]
        r['quality_string'] = line[10]
        r['mapq'] = int(line[4])

        # maybe add in other fields later...this should be enough to filter for now...
        
        r.append()

    def filter_reads_by_column(self,flt_function):
        # iterate over rows with an arbitrary function (which takes the row object as an argument); delete any rows for which function returns 'True'
        for r in range(0,len(read_table)):
            if (flt_function(read_table[r])):
                read_table.removeRows[r]

    def filter_unique_reads(self):
        pass
    

class GenomeTrack():
    def __init__(self,parent_table,species=None,genome_fasta=None,data_file=None,data_format=None,tag_limit=None,name=None):
        self.parent_table = parent_table
        self.species = species
        self.genome_fasta = genome_fasta
        self.data_file = data_file
        self.data_format = data_format
        self.tag_limit = tag_limit
        self.normalized = False
        self.array_tbl_idx = len(self.parent_table.listArrays)
        if not name:
            self.name = "array_" + self.array_tbl_idx
        else:
            self.name = name
        self.group = self.parent_table.hf5_object.createGroup(self.parent_table.array_group,self.name,'Tag Count Arrays ' + self.name)
        self.chromosomes = {}

        # build 0-filled arrays of genome
        if not ((species or genome_fasta) or group):
            raise TypeError("GenomeTrack() requires either a species specification or a genome")
        elif self.species:
            self.chromosomes = self.loadSpeciesArrays()
        elif self.genome_fasta:
            self.chromosomes = self.buildArraysFromFasta()

        # Populate Arrays
        if data_file and data_format.lower() == "wig":
            self.parse_wig(data_file)
        elif data_file and data_format.lower() == "bedgraph":
            self.parse_bedgraph(data_file)

    def add_tag(self,chromosome,loc,strand):
        if (not self.tag_limit) or (self.chromosomes[chromosome][strand][loc] < self.tag_limit):
            self.chromosomes[chromosome][strand][loc] += 1

    def return_extended(self,width):
        ext_chrs = {}
        for (chr,a_obj) in self.chromosomes.items():
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
            ext_chrs[chr] = ext_a_c_m + ext_a_c_p
            print "done extending - %f elapsed" % (time() - c_st)

    def reset_norm(self):
        self.norm = False

    def count_normalize(self,nfactor):
        if self.norm:
            print "WARNING: already normalized; use reset_norm() to dismiss"
        n_mult = nfactor / get_total_positions()
        for (chromosome,a_obj) in self.chromosomes.items():
            chromsome *= n_mult

    def window_normalize(self,window_size,function):
        pass
        
    def parse_wig(self,data_file):
        pass

    def parse_bedgraph(self,data_file):
        pass

    def outputBedGraph(self,outfile='out.bdg',step='10',description="User-Generated Track",viewmax=200):

        # TEST
        
        definition_line = "track type=bedGraph name=%s description=%s\
        visibility=full color=0,0,255 altColor=255,0,0\
        priority=priority autoScale=off alwaysZero=on\
        gridDefault=off maxHeightPixels=128:128:128\
        graphType=bar viewLimits=0:%d\
        yLineMark=0.0 yLineOnOff=on\
        windowingFunction=maximum smoothingWindow=off\
        " % (self.data_file,description,viewmax)
        out = open(outfile,'w')
        for (chromosome,vals) in self.chromosomes.items():
            print >> sys.stderr, "Writing Chr %s" % (chromosome,)
            for x in range(0,len(vals),step):
                total = 0
                if not (x+step >= len(vals)):
                    total = sum(vals[out_idx:out_idx+step])
                    print >> outfile, "%s\t%d\t%s" % (chromosome,x,total)
                else:
                    total = sum(vals[out_idx:len(vals) - 1])
                    print >> outfile, "%s\t%d\t%s" % (chromosome,x,total)
                    break
        outfile.close()
            
    def loadSpeciesArrays(self):
    # Get known species fasta from SPECIES_FASTAS dict, use to run buildArraysFromFasta()
        self.genome_fasta = SPECIES_FASTAS[self.species] 
        return self.buildArraysFromFasta()

    def buildArraysFromFasta(self):
    # Create empty arrays, one per chromosome.  Get length information from FastaDB.
        fa = FastaDB()
        fa.openFastaFile(self.genome_fasta)
        chr_name_size = max(map(len,fa.keys()))
        chr_arrays = {}
        for (chr,s) in fa.items():
            chr_arrays[chr] = {}
            chr_arrays[chr]["+"] = np.zeros(int(s['length']),dtype=np.uint32)
            chr_arrays[chr]["-"] = np.zeros(int(s['length']),dtype=np.uint32)
            total_positions += int(s['length'])

        for chr,a in chr_arrays.items():
            chr_carrays[chr] = {"-":None,"+":None}
            
            chr_carrays[chr]["+"] = self.parent_table.hf5_object.createCArray(self.group,"chr"+chr+"_plus",UInt32Atom(),(len(a['+']),),title=chr + "Plus Strand",filters=FILTERS)
            chr_carrays[chr]["-"] = self.parent_table.hf5_object.createCArray(self.group,"chr"+chr+"_minus",UInt32Atom(),(len(a['-']),),title=chr + "Minus Strand",filters=FILTERS)

            chr_carrays[chr]["+"][0:] = a['+']
            chr_carrays[chr]["-"][0:] = a['-']

        return chr_arrays

    def get_total_positions(self):
        tot = 0
        for arr in self.chromosomes.values():
            tot += len(arr)
        return tot
        

    def __get_item__(self,name):
        return self.chromsomes[name]
        
class chipExpt():
    #
    # Front-end class for interaction with ChIP data
    # 
    def __init__(self,group,name="NoName",chrs=DMEL_CHRS_EU,chr_ext='ext'):
        self.chr_dict = {}
        self.group = group
        self.name = name
        arrays = group._v_children.items()
        if chrs:
            for ch in chrs:
                chr_found = 0
                for (name,node) in arrays:
                    if re.search(chr_ext + "\w*_(chr)?"+ch+"$",name):
                        #print (ch,name,node.name)
                        chr_found = 1
                        self.chr_dict[ch] = node
                        break
                if chr_found == 0:
                    print "no chr found for %s" % (ch,)
        else:
            for n in group._f_walkNodes():
                ch_s = re.search("chr(.+)",n.name)
                self.chr_dict[ch_s.group(1)] = n
                break
                    
    def __getitem__(self,item):
        return self.chr_dict[item]

    def write_wig(self,chr_convert=None,span=10,outfile=None):
        if outfile == None:
            outfile = self.name + ".wig"
        wigout = open(outfile,"w")
        print >> wigout, "track type=wiggle_0 name=%s description=%s visibility=full autoScale=off maxHeightPixels=100:50:20" % (self.name,self.name)
        for (chrom,chr_arr) in self.chr_dict.items():
            print "Writing %s..." % (chrom,),
            name = chrom
            if chr_convert != None:
                name = chr_convert[chrom]
            print >> wigout, "variableStep chrom=%s span=%d" % (name,span)
            for i in np.arange(1,len(chr_arr),step=span):
                print >> wigout, "%d\t%f" % (i,np.sum(chr_arr[i:i+span])/float(span))
            print "done!"
        wigout.close()
