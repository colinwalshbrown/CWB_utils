#!/usr/bin/env python
"""
Module for handling sequence alignments.  
"""

import re
import os
import sys
import sqlite3
import pickle
import fasta_subseq_2 as fs2
import tables as tb
import numpy as np
import subprocess as sub
import tempfile as tmp
import seq_plotmethods as sp

BLAST_EXE = "/Users/cwbrown/src/blast3/blastn"
DMEL_GENOME = "/Users/cwbrown/Data/genomes/dmel.fa"
DB_DIR = os.environ["HOME"] + "/scripts/python_modules/alignment_DBs/"
ALN_H5 = os.environ["HOME"] + "/scripts/python_modules/aln_h5s/"
GENOME_DIR = os.environ["HOME"] + "/Data/genomes/"
SPP_FASTAS = {'dmel':GENOME_DIR + "dmel.fa",
              'dsim':GENOME_DIR + "dsim.fa",
              'dsec':GENOME_DIR + "dsec.fa",
              'dyak':GENOME_DIR + "dyak.fa",
              'dere':GENOME_DIR + "dere.fa",
              'dana':GENOME_DIR + "dana.fa",
              'dpse':GENOME_DIR + "dpse.fa",
              'dper':GENOME_DIR + "dper.fa",
              'dwil':GENOME_DIR + "dwil.fa",
              'dvir':GENOME_DIR + "dvir.fa",
              'dmoj':GENOME_DIR + "dmoj.fa",
              'dgri':GENOME_DIR + "dgri.fa",
              'dmir':GENOME_DIR + "dmir.fa",
              'dsub':GENOME_DIR + "dsub.fa",
              'dlow':GENOME_DIR + "dlow.fa",
              'daff':GENOME_DIR + "daff.fa",
              'dath':GENOME_DIR + "dath.fa",
              'dalg':GENOME_DIR + "dalg.fa"}
              

class Align:
    
    def __init__(self,aln_file=None,aln_format="fasta",aln_string=None):
        self.length = 0
        self.alignSeqs={}
        if aln_file:
            if aln_format == 'fasta':
                self.parseFasta(aln_file)
        if aln_string:
            self.parseString(aln_string)

    def parseString(self,string):
        name = ""
        seq = ""
        description = ""
        #print string.split("\n")
        for ln in string.split('\n'):
            if len(ln) < 1:
                continue
            if ln[0] == ">":
                if (name and seq):
                    self.addSeq(name,seq,description=description)
                n_search = re.search(">(\S+)\s*(.*)$",ln)
                name = n_search.group(1)
                description = n_search.group(2)
                seq = ""
            else:
                seq += ln
        if (name and seq):
            self.addSeq(name,seq,description=description)

    def parseFasta(self,fasta_file):
        name = ""
        seq = ""
        for ln in open(fasta_file):
            if ln[0] == ">":
                if (name and seq):
                    self.addSeq(name,seq)
                name = ln[1:-1]
                seq = ""
            else:
                seq += ln[:-1]
        if (name and seq):
            self.addSeq(name,seq)


    def getSeqDict(self):
        return self.alignSeqs

    def addSeq(self,name,seq,parent=None,parentStart=None,parentEnd=None,parentStrand=None,description=None):
        newSeq = AlignedSeq(seq,parent=parent,parentStart=parentStart,parentEnd=parentEnd,parentStrand=parentStrand,description=description)
        self.alignSeqs[name] = newSeq;
        if(self.length == 0):
            self.length = len(newSeq)
        elif(self.length != len(newSeq)):
            raise IOError("sequence " + name + " is not the correct length")
        
    def AlnSlice(self,start,end,seqs=[]):
        """
        create a new Align object from slices of sequences specificied by argument seq
        """
        if len(seqs) == 0:
            seqs = self.alignSeqs.keys()
        
        sliceAln = Align()
        for seqName in seqs:
            fullseq = self.alignSeqs[seqName]
            sliceSeq = fullseq[start:end]
            parentStart = None
            parentEnd = None
            parentStrand = None   
            if fullseq.parentStart:
                if fullseq.parentStrand == "-":
                    parentStart = fullseq.getParentCoord(end)
                else:
                    parentStart = fullseq.getParentCoord(start)
            if fullseq.parentEnd:
                if fullseq.parentStrand == "-":
                    parentEnd = fullseq.getParentCoord(start)
                else:
                    parentEnd = fullseq.getParentCoord(end)
            sliceAln.addSeq(seqName,sliceSeq,parent=fullseq.parent,parentStart=parentStart,parentEnd=parentEnd,parentStrand=fullseq.parentStrand)
        return sliceAln

    def AlnSliceFromUnaligned(self,seqName,start,end,seqs=[]):
        """
        creates an alignment slice from Align using start and end specified with unaligned
        coordinates from sequence seqName.  NOTE: includes ALL gapped columns from 
        alignment in between specified endpoints, even if they are gapped in the aligned
        sequence of seqName
        """
        if len(seqs) == 0:
            seqs = self.alignSeqs.keys()
   
        startAln = self.alignSeqs[seqName].getMapCoord(start)
        endAln = self.alignSeqs[seqName].getMapCoord(end)
        return self.AlnSlice(startAln,endAln,seqs)

    def AlnSliceFromAlnColumns(self,seqName,start,end,seqs=[]):
        """
        make alignment slice from seqName unaligned coordinates start and end; ONLY 
        include columns which don't have gaps in seqName 
        """
        if len(seqs) == 0:
            seqs = self.alignSeqs.keys()

        sliceAln = Align()
        sliceSeqs = dict.fromkeys(seqs,"")
        map = self.alignSeqs[seqName].getMap()
        for x in range(start,end):
            mapCoord = map[x]
            for seqName in seqs:
                sliceSeqs[seqName] += self.alignSeqs[seqName][mapCoord]
        for seqName in seqs:
            sliceAln.addSeq(seqName,sliceSeqs[seqName])
        return sliceAln

    def GetMatchLine(self,char="*"):
        matches = ""
        for x in range(0,self.length):
            curChar = ""
            for seq in self.alignSeqs.values():
                if (curChar == "") and (seq[x] != curChar):
                    curChar = seq[x].upper()
                elif (seq[x].upper() != curChar):
                    matches += " "
                    curChar = ""
                    break
            if curChar != "":
                matches += char
        return matches

    def GetNoGapLine(self,gapChar="-",char=".",maxGaps=0):
        matches = ""
        match_cnt = 0
        for x in range(0,self.length):
            gaps = 0
            for seq in self.alignSeqs.values():
                #if (curChar == "") and (seq[x] != curChar):
                #    curChar = seq[x].upper()
                if (seq[x].upper() == gapChar):
                    gaps += 1
            if gaps <= maxGaps:
                matches += char
                match_cnt += 1
            else:
                matches += " "
        print "matching positions: %d" % (match_cnt,)
        return matches

    
    def PrintClustal(self,lineLength,outHandle=sys.stdout):
        alnPtr = 0
        nameLim = 20
        while alnPtr < self.length:
            end = 0
            slice = None
            if (alnPtr + lineLength > self.length):
                slice = self.AlnSlice(alnPtr,alnPtr + (self.length - alnPtr))
            else:
                slice = self.AlnSlice(alnPtr,alnPtr + lineLength)
            for (seqName,seq) in slice.alignSeqs.iteritems():
                name = ""
                if len(seqName) > nameLim:
                    name = seqName[0:nameLim]
                else:
                    name = seqName + ((nameLim - len(seqName)) * " ")
                print >> outHandle, "%s\t%s" % (name,str(seq))
            print >> outHandle, "%s\t%s" % ((" " * nameLim),slice.GetMatchLine())
            alnPtr += lineLength

    def revcomAln(self):
        newAln = Align()
        newAln.length = self.length
        for (name,seq) in self.alignSeqs.items():
            newAln.alignSeqs[name] = seq.revcom()
        return newAln
            
    def transformSpCoord(self,sp1,sp2,coordinate):
        """
        translate parent coordinate from species 1 to species 2
        """
        sp1_aln_coord = self[sp1].getAlnCoordFromParent(coordinate)
        sp2_parent_coord = self[sp2].getParentCoord(sp1_aln_coord)
        return sp2_parent_coord

    def getSegregatingSites(self,include_gaps=False):
        sites = []
        seq_names = self.alignSeqs.keys()
        for i in range(self.length):
            site_bases = []
            for n in seq_names:
                site_bases.append(self.alignSeqs[n][i])
            if not include_gaps:
                site_bases = [j for j in site_bases if j != '-']
            if len(set(site_bases)) > 1:
                sites.append(i)
        return sites

    def getPrivateAlleles(self,seq,snps_only=False,filter_Ns=True):
        sites = []
        seq_names = self.alignSeqs.keys()
        for i in range(self.length):
            site_bases = []
            for n in seq_names:
                site_bases.append(self.alignSeqs[n][i])
            site_bases_notarget = [j for j in site_bases if j != self.alignSeqs[seq][i]]
            if snps_only and (("-" in site_bases) or ("N" in site_bases)):
                continue
            if filter_Ns and self.alignSeqs[seq][i] == "N":
                continue
            if len(site_bases_notarget) == (len(site_bases) - 1):
                sites.append(i)
        return sites


    def __str__(self):
        str_out = ""
        for seqName in self.alignSeqs:
            s = self.alignSeqs[seqName]
            if s.parent:
                str_out += ">%s\t%s:%d-%d:%s\n%s\n" % (seqName,s.parent,s.parentStart,s.parentEnd,s.parentStrand,str(s))
            else:
                str_out += ">%s\n%s\n" % (seqName,str(s))     
        return str_out

    def __getitem__(self,key):
        return self.alignSeqs[key]

    def __len__(self):
        return self.length

class AlignedSeq:    
    def __init__(self,seq,parent=None,parentStart=None,parentEnd=None,parentStrand=None,gap_char="-",description=None):
        self.seq = seq
        self.parent = parent
        self.parentStart = parentStart
        self.parentEnd = parentEnd
        self.parentStrand = parentStrand
        self.gap_char = gap_char
        self.description = description
        self.seqMap = self._makeSeqMap(seq,gap_char)
        
    def _makeSeqMap(self,aln_seq,gap_char):
        """
        build a mapping of coordinates in unaligned seq -> coordinates in aligned seq
        """
        map = []
        aln_count = 0
        for c in aln_seq:
            if c != gap_char:
                map.append(aln_count)
            aln_count += 1
        return map

    def getSliceFromUnaligned(self,start,end):
        """
        generate string for UNaligned sequence slice
        """
        slice = ""
        for x in range(start,end):
            slice += self.seq[self.seqMap[x]]
        return slice
    
    def getUnaligned(self):
        """
        generate entire unaligned sequence
        """
        unaligned = ""
        for x in self.seqMap:
            unaligned += self.seq[x]
        return unaligned

    def getLenUnaligned(self):
        return len(self.seqMap)
    
    def getMap(self):
        return self.seqMap

    def getDescription(self):
        return self.description
    
    def getMapCoord(self,coordinate):
        return self.seqMap[coordinate]

    def getUnalignedCoord(self,aln_coordinate):
        u_coord = None
        try:
            u_coord = self.seqMap.index(aln_coordinate)
        except ValueError:
            u_coord = self.seqMap.index(minGTZero(map(lambda x: x - aln_coordinate,self.seqMap)) + aln_coordinate)
            #if (u_coord - aln_coordinate) < 0:
            #    u_coord = len(self.seqMap) - 1
        return u_coord

    def getParentCoord(self,aln_coordinate):
        p_coord = None
        if self.parentStrand == "-":
            p_coord = self.parentStart + (self.getLenUnaligned() - self.getUnalignedCoord(aln_coordinate)) 
        else:
            p_coord = self.parentStart + self.getUnalignedCoord(aln_coordinate) 
        return p_coord

    def getParentCoordsMap(self,ret = 'array'):
        # return an array (or 'list') of map coordinates relative to starting position of the parent (e.g. chromosomal coordinates)
        map_arr = np.arange(0,self.getLenUnaligned(),dtype=np.uint64)
        p_coords = None
        if self.parentStrand == "-":
            ualn = np.repeat(np.array([self.getLenUnaligned()],dtype=np.uint64),self.getLenUnaligned())
            p_coords = (ualn - map_arr) + self.parentStart
        else:
            p_coords = map_arr + self.parentStart
        if ret == 'list':
            p_coords = list(p_coords)
        return p_coords

    def getAlnCoordFromParent(self,parentCoord):
        ualn_coord = None
        #print parentCoord
        #print self.parentStart
        #print self.parentEnd
        if parentCoord >= self.parentEnd:
            raise IndexError()
            parentCoord = self.parentEnd -1
        elif parentCoord <= self.parentStart:
            raise IndexError()
            parentCoord = self.parentStart + 1

        if self.parentStrand == "-":
            ualn_coord = self.getLenUnaligned() - (parentCoord - self.parentStart)
        else:
            ualn_coord = parentCoord - self.parentStart

        #print ualn_coord
        aln_coord = self.getMapCoord(ualn_coord) 
        return aln_coord

    def revcom(self):
        newStrand = None
        if self.parentStrand == '-':
            newStrand = "+"
        elif self.parentStrand == "+":
            newStrand = "-"
        newseq = AlignedSeq(revcomp(self.seq),self.parent,self.parentStart,self.parentEnd,newStrand)
        return newseq

    def blast_seq(self,genome=DMEL_GENOME,name="query_seq",top=5):
        all_hits = []
        fa_seq = ">%s\n%s" % (n,self.getUnaligned())
        temp_q = tmp.NamedTemporaryFile(suffix=".fa")
        temp_q_name = temp_q.name
        print >> temp_q, fa_seq
        temp_q.flush()
        blast_proc = sub.Popen([BLAST_EXE,genome,temp_q_name,"-mformat","3"],stdout=sub.PIPE,stderr=sub.PIPE)
        (out,err) = blast_proc.communicate()
        temp_q.close()
        #print out
        #print err
        hits = [x for x in out.split("\n") if (len(x) > 1) and (x[0] != "#")]
        comments = [x for x in out.split("\n") if (len(x) > 1) and (x[0] == "#")]
        for bl_ln in hits[:top]:
            bl_ln_sp = bl_ln.strip().split()
            bl = {'chr':bl_ln_sp[1],
                  'start':int(bl_ln_sp[20]),
                  'end':int(bl_ln_sp[21]),
                  'length':int(bl_ln_sp[21]) - int(bl_ln_sp[20]),
                  'strand':"+",
                  'name':name}
            if bl_ln[16] == "-1":
                bl['strand'] = "-"
            all_hits.append(bl)
        return all_hits

    def __getitem__(self,key):
        return self.seq[key]
 
    def __getslice__(self,start,end):
        return self.seq[start:end]
  
    def __len__(self):
        return len(self.seq)

    def __str__(self):
        return self.seq


class WholeGenomeAlign:
    
    def __init__(self,alnFile=None,alnh5DB=None,format=None,reindex=False,backend='h5'):
        self.alnFile = alnFile
        self.format = format
        self.alnDBs = None
        self.alnDBs = pickle.load(open(os.environ["HOME"] + "/scripts/python_modules/ALIGNMENT_INDEX_DBS"))
        self.alnh5DB = alnh5DB
        self.aln = None
        if (self.alnFile):
            print self.alnDBs
            if (alnFile in self.alnDBs.keys()) and not reindex:
                self.aln = WGAlignDB(self.alnDBs[self.alnFile])
            elif (alnFile in self.alnDBs.keys()) and reindex:
                print self.alnDBs[alnFile]
                try:
                    os.unlink(self.alnDBs[alnFile])
                except OSError:
                    print "Couldn't delete %s - file may be missing, will attempt overwrite" % self.alnDBs[alnFile]
                del self.alnDBs[alnFile]
                self.aln = self.parseWGAlign(self.alnFile,self.format)
                self.alnDBs[self.alnFile] = self.aln.DBname
                pickle.dump(self.alnDBs,open(os.environ["HOME"] + "/scripts/python_modules/ALIGNMENT_INDEX_DBS","w"))
            else:
                self.aln = self.parseWGAlign(self.alnFile,self.format)
                self.alnDBs[self.alnFile] = self.aln.DBname
                pickle.dump(self.alnDBs,open(os.environ["HOME"] + "/scripts/python_modules/ALIGNMENT_INDEX_DBS","w"))
        elif (self.alnh5DB):
            self.aln = WGAlignDB_h5(self.alnh5DB)
        

    def parseWGAlign(self,file,format):
        alnObj = None
        if format.lower() == "maf":
            alnObj = self.parseMAF(file)
        elif format.lower() == "fsa":
            alnObj = self.parseFSADir(file)
        return alnObj

    def parseFSADir(self,map_file):
        """
        parse a directory of FSA fasta (mfa) alignments; assumes existence of a map and genomes file in same dir as alignments
        """
        base_dir = os.path.dirname(map_file)
        DB_name = os.path.basename(base_dir)+"_WGAlign_DB"
        alnDB = WGAlignDB(DB_DIR + DB_name)
        genome_file = open(base_dir + '/' + 'genomes')
        genomes = genome_file.readline()[:-1].split()
        for map_line in open(map_file):
            map_spl = map_line[:-1].split()
            aln_file = map_spl[0]
            coords = map_spl[1:]
            s = 0
            aln_in = AlignIO(base_dir + '/' + aln_file + ".mfa")
            alnObj = aln_in.getAlign()
            for (name,aln) in alnObj.getSeqDict().items():
                seqname = name.split("_")[0]
                alnObj.alignSeqs[seqname] = aln
                alnObj.alignSeqs.pop(name)
            while ((s + 1)*3 <= len(coords)):
                seqname = genomes[s].split("_")[0]
                (chr,start,end,strand) = coords[4*s:4*s+4]
                print (chr,start,end,strand)
                alnObj.alignSeqs[seqname].parent = chr
                alnObj.alignSeqs[seqname].parentStart = int(start)
                alnObj.alignSeqs[seqname].parentEnd = int(end)
                alnObj.alignSeqs[seqname].parentStrand = strand
                s += 1
            alnDB.addAln(alnObj,None)
            print >> sys.stderr, "Added %s" % map_line
        alnDB.commitAdds()
        return alnDB

    def parseMAF(self,file):
        DB_name = ".".join(file.split("/")[-1].split(".")[:-1]) + "_WGAlign_DB"
        #print DB_DIR + DB_name
        alnDB = WGAlignDB(DB_DIR + DB_name)
        alnObj = None
        score = None
        s = 0
        for x in open(file):
            #print x
            # parse a MAF line into an Align object
            comm_re = re.search("^#",x)
            score_re = re.search("^a\s+score=(\S+)",x)
            seq_re = re.search("^s\s+(\S+)\.(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)",x)
            if comm_re:
                continue
            elif score_re:
                #print score_re.groups()
                score = float(score_re.group(1))
            elif seq_re and not alnObj:
                alnObj = Align()
                #print seq_re.groups()
                alnObj.addSeq(seq_re.group(1),seq_re.group(7),parent=seq_re.group(2),
                              parentStart=int(seq_re.group(3)),parentEnd=(int(seq_re.group(3))+int(seq_re.group(4))),
                              parentStrand=seq_re.group(5))
            elif seq_re and alnObj:
                #print seq_re.groups()
                alnObj.addSeq(seq_re.group(1),seq_re.group(7),parent=seq_re.group(2),
                              parentStart=int(seq_re.group(3)),parentEnd=(int(seq_re.group(3))+int(seq_re.group(4))),
                              parentStrand=seq_re.group(5))
            else:
                # add the object to the alignment dict
                s += 1
                alnDB.addAln(alnObj,score)
                alnObj = None
                score = None
            print s
            
        alnDB.commitAdds()
        return alnDB

    def getAllAlns(self):
        print self.aln
        alns = self.aln.DBcurs.execute("SELECT * FROM subaln")
        for a in alns:
            print pickle.loads(a["alnobj"].encode('utf-8'))

    # need to move this to be DB-type-specific (i.e. move this one to the WGAlignDB_sql class)

def WGAlignDB(filename,backend='h5',reindex=False):
        if backend == 'sql':
            return WGAlignDB_sql(filename)
        elif backend == 'h5':
            return WGAlignDB_h5(filename+".h5",reindex=reindex)

class WGAlignDB_h5:
    #
    # re-implement WGAlign using pytables for better retrieval speed
    #
    def __init__(self,h5,reindex=False):
        self.DBname = h5
        self.reindex = reindex
        self.reindex_sp = []
        if self.reindex:
            try:
                os.remove(h5)
            except OSError:
                print "file %s does not exist - creating" % (h5,)
        self.aln_table = None
        self.species_chrs = {}
        self.chr_key_arrays = {}
        self.built_chr_tabs = False
        self.species = []
        print self.DBname
        self.h5 = None
        if not self.reindex:
            self.h5 = tb.openFile(self.DBname)
            self._fillWGADBFromFile()
        else:
            self.h5 = tb.openFile(self.DBname,"w")

    def __del__(self):
        self.h5.close()

    def _fillWGADBFromFile(self):
        #print self.h5.root
        #print self.h5.root._v_children.items()
        self.aln_table = self.h5.root.aln_table
        for sp in self.h5.root._v_groups.keys():
            self.species_chrs[sp] = {}
            self.species.append(sp)
            sp_fasta = fs2.FastaDB()
            sp_fasta.openFastaFile(SPP_FASTAS[sp])
            chrs = sorted(sp_fasta.keys())#.sorted()
            self.chr_key_arrays[sp] = [None] * len(chrs)
            for chrom in self.h5.root._v_groups[sp]._v_children.keys():
                ch = chrom.replace("chr","")
                #print (ch,chrom)
                self.species_chrs[sp][ch] = self.h5.root._v_groups[sp]._v_children[chrom]
            for (i,chrom) in enumerate(chrs):
                #print self.species_chrs
                try:
                    self.chr_key_arrays[sp][i] = self.species_chrs[sp][chrom]
                except:
                    print >> sys.stderr, "WARNING: chromsome %s not found" % (chrom,)
        self.built_chr_tabs = True
        
    def addAln(self,subAlnObj,alnScore):
        if (not self.built_chr_tabs):
             self._buildAlnTables(subAlnObj,alnScore)
        #insert_len = len(subAlnObj.getSeqDict())
        insert = np.ndarray(shape=(subAlnObj.length,1),dtype=self.aln_tbl_dtype)
        for (sp,seq) in subAlnObj.getSeqDict().items():
            print "adding %d bases from %s" % (len(seq),sp)
            seq_ar = np.matrix(list(str(seq)),dtype=np.dtype('a1')).T
            print "built seq mtx"
            chr_ar = np.matrix([self.chr_key_arrays[sp].index([seq.parent,self.species_chrs[sp][seq.parent]])] * len(seq),dtype=np.dtype('u4')).T
            print "built chr mtx"
            par_coords = seq.getParentCoordsMap()
            #print par_coords[:1000]
            map_ar = np.matrix(np.zeros(len(seq_ar)),dtype=np.dtype('u8')).T
            print "par_coords length: %d" % len(par_coords)
            par_ix = 0
            ingap = False
            gap_st = [0,0]
            gap_end = [0,0]
            gap_len = 0
            par_gap_start = 0
            par_gap_end = 0
            # redo this in cython maybe? Kinda slow...
            for (i,b) in enumerate(seq_ar):
                if (b == seq.gap_char and not ingap):
                    par_gap_start = par_coords[par_ix - 1]
                    gap_st[0] = i
                    gap_len = 1
                    ingap = True
                elif (b == seq.gap_char and ingap):
                    gap_len += 1
                elif ingap:
                    if (gap_len) <= 1:
                        map_ar[gap_st[0]] = par_gap_start
                        
                    gap_st[1] = (gap_st[0] + int(gap_len / 2))
                    gap_end[0] = gap_st[1]
                    gap_end[1] = i
                    par_gap_end = par_coords[par_ix]
                    map_ar[gap_st[0]:gap_st[1]] = np.repeat(np.matrix([par_gap_start],dtype=np.dtype('u8')),(gap_st[1] - gap_st[0])).T
                    map_ar[gap_end[0]:gap_end[1]] = np.repeat(np.matrix([par_gap_end],dtype=np.dtype('u8')),(gap_end[1] - gap_end[0])).T
                    map_ar[i] = par_gap_end
                    par_ix += 1
                    ingap = False
                else:
                    map_ar[i] = par_coords[par_ix]
                    par_ix += 1

            if (ingap):
                gap_st[1] = len(seq_ar) 
                map_ar[gap_st[0]:] = np.repeat(np.matrix([par_gap_start],dtype=np.dtype('u8')),(gap_st[1] - gap_st[0])).T
            
            print "built map mtx"
            sp_idx = self.species.index(sp)
            insert [sp]['base'] = seq_ar
            print "added base to insert table"
            insert [sp]['chr_key'] = chr_ar
            print "added chr_key to insert table"
            insert [sp]['position'] = map_ar
            print "added position to insert table"
        nxti = len(self.aln_table)
        newrs = nxti + len(insert)
        print "updating chr maps..."
        self._update_chrs(insert,nxti,newrs)
        print "adding remaining rows..."
        self.aln_table.append(insert)
        print "done!"

    def _update_chrs(self,insert,nxti,newrs):
        for (i,r) in enumerate(range(nxti,newrs)):
            irow = insert[i]
            for s in self.species:
                base = self.chr_key_arrays[s][irow[s]['chr_key']][1][int(irow[s]['position'])]['base']
                self.chr_key_arrays[s][irow[s]['chr_key']][1][int(irow[s]['position'])] = [base,r]
  
    def _buildAlnTables(self,subAlnObj,alnScore):
        species = subAlnObj.getSeqDict().keys()
        dt_init = [] # description for alignment table - build as we read spp
        exp_size = 0
        self.species = species
        for (j,sp) in enumerate(species):
            # open fasta to get unaligned chromosome lengths & seqs
            sp_fasta = fs2.FastaDB()
            sp_fasta.openFastaFile(SPP_FASTAS[sp])
            chrs = sorted(sp_fasta.keys())
            # group for each species' chromosome arrays
            sp_grp = self.h5.createGroup(self.h5.root,sp,"%s Chromosomes" % (sp,))
            self.species_chrs[sp] = {}
            # use 1 byte index as chr identifier
            self.chr_key_arrays[sp] = [None] * len(chrs)
            for (i,ch) in enumerate(chrs):
                # build 2 x len(chr) array of (base,aligned_coord) pairs
                # print sp_fasta[ch]
                bases = np.matrix(list(sp_fasta[ch].getFullSeq()),dtype=np.dtype("a1")).T
                maps = np.matrix(np.zeros(len(sp_fasta[ch]) - 1),dtype=np.dtype("u8")).T
                flat_chr_arr = np.ndarray(shape=(len(maps),1),dtype=np.dtype([('base','a1'),('aln_map','u8')]))
                flat_chr_arr['base'] = bases
                flat_chr_arr['aln_map'] = maps 
                self.species_chrs[sp][ch] = self.h5.createTable(sp_grp,"chr" + ch,np.dtype([('base','a1'),('aln_map','u8')]))
                self.species_chrs[sp][ch].append(flat_chr_arr)
                self.species_chrs[sp][ch].flush()
                self.chr_key_arrays[sp][i] = [ch,self.species_chrs[sp][ch]]
                print "%s %s length = %d added to %s" % (sp,ch,len(bases),sp_grp)
                exp_size += len(sp_fasta[ch])
            # add a column to the align table description
            dt_init.append((sp,[('base','a1'),('chr_key','u4'),('position','u8')]))
        self.aln_tbl_dtype = np.dtype(dt_init)
        self.aln_table = self.h5.createTable(self.h5.root,'aln_table',self.aln_tbl_dtype,expectedrows=exp_size) # make the alignment table
        self.aln_table.flush()
        self.built_chr_tabs = True

    def getSliceSpeciesCoord(self,refsp,chrom,start,end,contiguous=True):
        chr_pos_base = self.species_chrs[refsp][chrom][start:end]
        chr_pos = chr_pos_base['aln_map']
        chr_base = chr_pos_base['base']
        ref_strand = None
        other_spp = list(self.species)
        other_spp.remove(refsp)
        aln_obj = Align()
        aln_objs = []
        aln_slice = None
        #print >> sys.stderr, chr_pos[0]
        if ((chr_pos[0] == 0 or chr_pos[-1] == 0) or (abs(chr_pos[-1] - chr_pos[0]) > 20 * len(chr_base))): #preeeeeetty kludgy...could adjust
            raise IndexError("Chimeric sequence or out-of-bounds for %s %s %d->%d" % (refsp,chrom,start,end))
        elif chr_pos[-1] > chr_pos[0]:
            try:
                ref_fwd = True
                aln_slice = self.aln_table[chr_pos[0]:chr_pos[-1]]
            except ValueError:
                print >> sys.stderr, "Coordinates [%d:%d] failed for array [0:%d]" % (chr_pos[0],chr_pos[-1],len(self.aln_table))
                print >> sys.stderr, "%s %s %d->%d" % (refsp,chrom,start,end)
        else:
            ref_fwd = False
            aln_slice = self.aln_table[chr_pos[-1]:chr_pos[0]]
        aln_obj.addSeq(refsp,"".join(aln_slice[refsp]['base']),parent=chrom,parentStart=start,parentEnd=end,parentStrand="+")


        for sp in other_spp:
            sp_slice = aln_slice[sp]
            sp_seq = "".join(sp_slice['base'])

            # deal with aligned seqs that span multiple chromosomes
            sp_chr = sp_slice['chr_key']
            sp_chimeric = False
            #print "%s %d %d %d" % (sp,sp_chr[0],sp_chr[-1],len(sp_chr))
            if (sum(sp_chr - sp_chr[0]) != 0):
                print "WARNING: Chimeric Sequence for %s:%s:%d-%d" % (refsp,chrom,start,end)
                sp_chimeric = True
                
            sp_st = sp_slice['position'][0]
            sp_end = sp_slice['position'][-1]
            sp_strand = None
            if (sp_st >= sp_end):
                if ref_fwd:
                    sp_strand = "-"
                    hold = sp_st
                    sp_st = sp_end
                    sp_end = hold
                else:
                    sp_strand = "+"
                    hold = sp_st
                    sp_st = sp_end
                    sp_end = hold

            elif (sp_st < sp_end):
                if ref_fwd:
                    sp_strand = "+"
                else:
                    sp_strand = "-"
            parent_chr_srch = re.search("chr(.+)",self.chr_key_arrays[sp][sp_chr[0]].name)
            parent_chr = parent_chr_srch.group(1)
            aln_obj.addSeq(sp,sp_seq,parent=parent_chr,parentStart=sp_st,parentEnd=sp_end,parentStrand=sp_strand)
            if contiguous and sp_chimeric:
                raise IndexError("Chimeric sequence found for %s alignment to %s:%s:%d-%d" % (sp,refsp,chrom,start,end))

        aln_objs.append(aln_obj)
            
        return aln_objs

    def translateSpeciesCoords(self):
        pass
        
    def commitAdds(self):
        self.aln_table.flush()

    class ChimericAlignedSeq():
        def __init__(self):
            self.chrs = []
            self.seqMap = []
            self.seqs = {}

        def addSeq(self,seq,parent,parentStart,parentEnd,parentStrand):
            newseq = AlignedSeq(seq,parent,parentStart,parentEnd,parentStrand)
            self.seqs[parent] = newseq
            self.chrs.append(parent)
            self.seqMap.extend(zip([parent] * len(newseq.seqMap),newseq.seqMap))

        def __getslice__(self,start,end):
            sl_seq = ""
            sl = self.seqMap[start:end]
            for b in sl:
                sl_seq += self.seqs[b[0]][b[1]]
            
                
    

class WGAlignDB_sql:
    # Wrapper class for whole genome alignments - stored in sqlite DB
    #
    def __init__(self,DB,reindex=False):
        self.DBname = DB
        self.reindex = reindex
        self.reindex_sp = []
        self.DBconn = sqlite3.connect(DB)
        self.DBconn.row_factory = sqlite3.Row
        self.DBcurs = self.DBconn.cursor()
        # initialize subAln table
        
        self.DBcurs.execute("""DROP TABLE IF EXISTS subaln""")
        self.DBcurs.execute("""CREATE TABLE IF NOT EXISTS subaln (subaln_key INTEGER PRIMARY KEY AUTOINCREMENT,
                                                                  aln_score FLOAT,
                                                                  alnobj STR)""")
    def addAln(self,subAlnObj,alnScore):
        # import an alignment object
        self.DBcurs.execute("INSERT INTO subaln VALUES (?,?,?)",(None,alnScore,pickle.dumps(subAlnObj)))
        insertID = self.DBcurs.lastrowid
        seqDict = subAlnObj.getSeqDict()
        for sp in seqDict.keys():
            if (self.reindex) and (sp not in self.reindex_sp):
                self.DBcurs.execute("""DROP TABLE IF EXISTS %s""" % sp)
                self.reindex_sp.append(sp)
            self.DBcurs.execute("""CREATE TABLE IF NOT EXISTS %s (%s INTEGER PRIMARY KEY AUTOINCREMENT,
                                                                 chr TEXT,
                                                                 start INT,
                                                                 end INT,
                                                                 subaln_key INT)""" % (sp,sp + "_key"))
            self.DBcurs.execute("CREATE INDEX IF NOT EXISTS %s_idx ON %s (chr,start,end)" % (sp,sp))
            self.DBcurs.execute("INSERT INTO %s VALUES (?,?,?,?,?)" % sp,(None,seqDict[sp].parent,
                                                                          seqDict[sp].parentStart,
                                                                          seqDict[sp].parentEnd,
                                                                          insertID))
            #print ("INSERT INTO %s VALUES (-,%s,%s,%s,%s)",(sp,seqDict[sp].parent,
            #                                                seqDict[sp].parentStart,
            #                                                seqDict[sp].parentEnd,
            #                                                insertID))


    def commitAdds(self):
        self.DBconn.commit()

    def getSliceSpeciesCoord(self,species,chr,start,end,contiguous=False):

        alns = self.DBcurs.execute("""select alnobj from subaln where subaln_key in 
                                   (select subaln_key from %s where ((chr = ?) and
                                   (((start between ? and ?) or 
                                   (end between ? and ?)) or
                                   ((? between start and end) or
                                   (? between start and end)))))""" % species,(chr,start,end,start,end,start,end))
        alnObjs = map((lambda x: pickle.loads(x['alnobj'].encode('utf-8'))),alns.fetchall())
        sortedAlnObjs = sorted(alnObjs,cmp=lambda x,y: cmp(x[species].parentStart,y[species].parentStart))
        alnSlices = []

        for aln in sortedAlnObjs:

            slice_start = None
            slice_end = None
            startBtw = between(start,aln[species].parentStart,aln[species].parentEnd)
            endBtw = between(end,aln[species].parentStart,aln[species].parentEnd)

            if startBtw and endBtw:
                if aln[species].parentStrand == "-":
                    ualn_start = aln[species].getLenUnaligned() - (end - aln[species].parentStart)
                    ualn_end = aln[species].getLenUnaligned() - (start - aln[species].parentStart)
                else:
                    ualn_start = start - aln[species].parentStart
                    ualn_end = end - aln[species].parentStart
                slice_start = aln[species].getMapCoord(ualn_start)
                slice_end = aln[species].getMapCoord(ualn_end)
            elif startBtw and not endBtw:
                if aln[species].parentStrand == "-":
                    ualn_end = aln[species].getLenUnaligned() - (start - aln[species].parentStart)
                    slice_end = aln[species].getMapCoord(ualn_end)
                    slice_start = 0
                else:
                    ualn_start = start - aln[species].parentStart
                    slice_start = aln[species].getMapCoord(ualn_start)
                    slice_end = aln.length - 1
            elif not startBtw and endBtw:
                if aln[species].parentStrand == "-":
                    ualn_start = aln[species].getLenUnaligned() - (end - aln[species].parentStart)
                    slice_start = aln[species].getMapCoord(ualn_start)
                    slice_end = aln.length - 1
                else:
                    ualn_end = end - aln[species].parentStart
                    slice_end = aln[species].getMapCoord(ualn_end)
                    slice_start = 0
              
            if slice_start or slice_end:
                alnSlices.append(aln.AlnSlice(slice_start,slice_end))
            else:
                alnSlices.append(aln)
            
        return alnSlices

   
class AlignIO:
    
    def __init__(self,file,format="fasta"):
        self.file = open(file)
        self.format = format
        self.align = self._readAlign(self.file,self.format)

    def _readAlign(self,file,format):
        if (format == 'fasta'):
            return self._readFastaAlign(file)
        else:
            raise  IOError("Unknown format: " + format)

    def _readFastaAlign(self,file):
        newAlign = Align()
        curSeqName = ""
        curSeq = ""

        for line in file:
            line_match = re.match("^>\s*(\S+)(\s+(.*))?",line)
            if ((line_match) and (curSeqName == "")):
                curSeqName = line_match.group(1)
            elif (line_match):
                newAlign.addSeq(curSeqName,curSeq)
                curSeqName = line_match.group(1)
                curSeq = ""
            else:
                curSeq += line[:-1]
        newAlign.addSeq(curSeqName,curSeq)
        return newAlign
    
    def getAlign(self):
        return self.align

class WGAlign_CNDSLice:
    
    def __init__(self, st_dir):
        self.st_dir = st_dir
        self.mapfile = open(st_dir + "/map")
        self.genomesfile = open(st_dir + "/genomes")
        
        self.genomes = self.genomesfile.readline().strip().split()
        self.map = []
        # use the alignment's Mercator map file to build index of aligned regions for each genome
        for ml in self.mapfile:
            mlsp = ml.strip().split()
            mp = {}
            for (i,sp) in enumerate(self.genomes):
                sp_mp = {}
                if mlsp[4*i + 1] == "NA":
                    (sp_mp['chr'],sp_mp['start'],sp_mp['end'],sp_mp['strand']) = (None,-1,-1,None)
                else:
                    sp_mp['chr'] = mlsp[4*i + 1]
                    sp_mp['start'] = int(mlsp[4*i + 2])
                    sp_mp['end'] = int(mlsp[4*i + 3])
                    sp_mp['strand'] = mlsp[4*i + 4]
                mp[sp] = sp_mp
                
            self.map.append(mp)
        
    def getSliceSpeciesCoord(self,species,chrom,start,end,fill_missing=False):
        """
        Use Colin Dewey's sliceAlignment program to get slice of a whole-genome alignment given one species' coordinates.

        Return an array of CWB Align objects (targets spanning multiple alignment blocks are split into separate alignment objects)

        fill_missing = True -> create a sequence of all "-"s for species not present in alignment

                       False -> Raise an AlignError exception if run fails or no alignment found
        """
        t_dict = {'chr' : chrom, 'start' : start, 'end' : end}
        alns = []
        map_hits = []
        
        for (i,mp) in enumerate(self.map):
             # find alignment blocks overlapping target sequence
            if sp.overlap(t_dict,mp[species]) != (0,0):

                # trim target sequence to only overlap a single alignment block
                ol_st = max(t_dict['start'],mp[species]['start'])
                ol_end = min(t_dict['end'],mp[species]['end'])

                # Call  to get alignment slice for trimmed target
                runalnsl = sub.Popen(["sliceAlignment",self.st_dir,species,chrom,str(ol_st),str(ol_end),"+"],stdout=sub.PIPE,stderr=sub.PIPE)
                (aln_str,aln_err) = runalnsl.communicate()

                # Raise AlignError exception if search fails
                if (len(aln_str) < 1) or (re.search("Warning",aln_err)):
                    raise AlignError("sliceAlignment run error for %s %s:%d->%d: %s" % (species,chrom,ol_st,ol_end,aln_err))

                # Make an Align object for result
                aln = Align(aln_string=aln_str)

                # Fill in parent sequence information for aligned seqs
                for spp in aln.alignSeqs.keys():
                    seq = aln[spp]
                    d = seq.getDescription()
                    #print d
                    d_search = re.search("(\S+?)\:(\d+)\-(\d+)(\+|\-)",d)
                    if (d_search):
                        seq.parent=d_search.group(1)
                        seq.parentStart=int(d_search.group(2))
                        seq.parentEnd=int(d_search.group(3))
                        seq.parentStrand=d_search.group(4)
                    elif not fill_missing:
                        # raise exception if no sequence for current spp
                        raise AlignError("No sequence for %s in alignment %s %s:%d->%d" % (spp,species,chrom,start,end))
                alns.append(aln)
        if len(alns) < 1:
            # raise exception if no alignment blocks found
            raise AlignError("No alignment blocks found for %s %s:%d->%d" % (species,chrom,start,end))
        return alns
                

class AlignError(Exception):
    def __init__(self, value):
        self.val = value
    def __str__(self):
        return repr(self.val)

# package methods
#
#
def between(x,y,z):
        if (y <= x) and (z >= x):
            return True
        else:
            return False

def minGTZero(list):
    min = list[0]
    for x in list:
        if x < min and x >= 0:
            min = x
        elif min < 0 and x > 0:
            min = x
    return min

def revcomp(seq):
    """
    return reverse complement of a string
    """ 
    newseq = ""
    comp = {'A':'T',
            'T':'A',
            'G':'C',
            'C':'G',
            'a':'t',
            't':'a',
            'c':'g',
            'g':'c'}
    for i in range(0,len(seq)):
       try:
           newseq = comp[seq[i]] + newseq
       except KeyError:
           newseq = seq[i] + newseq

    return newseq
