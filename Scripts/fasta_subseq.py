#!/usr/bin/env python
"""
get subsequence of a fasta file. Create index for future use if necessary.
"""

import sys
import re

class FastaDB:

    def __init__(self):
        pass

    def makeIndex(self):

        self.index = {}
        curseq = "";
        offset = 0;

        for line in self.file:

            offset += len(line)

            line_match = re.match("^>(\S+)\s+(.+)",line)
            
            if(line_match):

                curseq = line_match.groups()[0]

                self.index[curseq] = {}

                for descriptor in line_match.groups()[1].split("; "):
                    desc = descriptor.split("=")
                    self.index[curseq][desc[0]] = desc[1]
    
                self.index[curseq]['sequence'] = FastaSeq(offset,self.file,int(self[curseq]['length']))

                curseq = line_match.groups()[0]

            else:
                self.index[curseq]['sequence'].setLineLen(len(line) - 1)
                

    def openFastaFile(self,file):
        
        self.file = open(file)

        self.makeIndex()

    def __getitem__(self,key):
        return self.index[key]

class FastaSeq():
    
    def __init__(self, offset, file, length):
        self.offset = offset
        self.file = file
        self.line_len = 0
        self.length = length

    def __getslice__(self,start,end):
        start -= 1
        nret_start = start / self.line_len
        nret_end = (end - start) / self.line_len
        
        self.file.seek(self.offset + start + nret_start)
        slice = self.file.read((end - start) + nret_end)
        
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
