#!/usr/bin/env python
"""
get subsequence of a fasta file. Create index for future use if necessary.
"""

import sys
import os
import re
import shelve


class FastaDB():

    def __init__(self,fasta_file=None,reindex=False):
        if (fasta_file):
            self.openFastaFile(fasta_file,reindex=reindex)
            
    def makeIndex(self):

        print "rebuilding index..."

        self.index = {}
        curseq = ""
        offset = 0
        length = 0

        for line in open(self.file):

            offset += len(line)
            
            line_match = re.match("^>\s*(\S+)(\s+(.+))?",line)
            
            if(line_match):
                if (curseq) and (self.index[curseq]['sequence'].length == -1):
                    self.index[curseq]['sequence'].length = length
                    length = 0
                curseq = line_match.groups()[0]

                self.index[curseq] = {}
                
                if line_match.groups()[1]:
                    for descriptor in line_match.groups()[1].split("; "):
                        desc = descriptor.split("=")
                        if len(desc) < 2:
                            continue
                        self.index[curseq][desc[0]] = desc[1]
                
                if self.index[curseq].has_key('length'):
                    self.index[curseq]['sequence'] = FastaSeq(offset,self.file,int(self.index[curseq]['length']),annotations=line[:-1])
                else:
                    self.index[curseq]['sequence'] = FastaSeq(offset,self.file,-1,annotations=line[:-1])
                
                curseq = line_match.groups()[0]

            else:
                self.index[curseq]['sequence'].setLineLen(len(line) - 1)
                length += len(line) - 1

    def openFastaFile(self,fasta_file,reindex=False):
        shelfdir = os.environ["HOME"] + "/scripts/python_modules/fasta_subseq_indices.shelf"
        idxShelf = shelve.open(shelfdir)

        if idxShelf.has_key(fasta_file) and not reindex: 
            self.index = idxShelf[fasta_file]
        else:
            self.file = fasta_file
            self.makeIndex()
            idxShelf[fasta_file] = self.index
        idxShelf.close()

    def __getitem__(self,key):
        return self.index[key]['sequence']

    def keys(self):
        return self.index.keys()

    def items(self):
        return self.index.items()
    
    def values(self):
        return self.index.values()

class FastaSeq():
    
    def __init__(self, offset, file, length, annotations=None):
        self.offset = offset
        self.file = file
        self.line_len = 0
        self.length = length
        self.annotations = annotations

    def getFullSeq(self):
        return self.__getslice__(0,self.length)

    def __getslice__(self,start,end):
        if self.length == -1:
            pass
        elif (start < 0) or (end > (self.length)):
            raise FastaSubseqError("index out of range")
        
        if not end:
            end = self.length
        
        if not start:
            start = 0
        else:
            start -= 1

        fileObj = open(self.file)
        nret_start = start / self.line_len
        nret_end = (end / self.line_len) - nret_start
        
        fileObj.seek(self.offset + start + nret_start)
        slice = fileObj.read((end - start) + nret_end)
        
        trim_slice = re.sub("\n","",slice)
        return trim_slice

    def __len__(self):
        return self.length

    def setLineLen(self,line_len):
        if(self.line_len == 0):
            self.line_len = line_len

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

class FastaSubseqError(Exception):
    def __init__(self, value):
        self.val = value
    def __str__(self):
        return repr(self.val)


