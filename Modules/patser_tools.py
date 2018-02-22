#!/usr/bin/env python
"""
patser_tools.py

Module for creating and manipulating annotations based from patser
"""

import re
import os
import cStringIO
import annotation_track
import subprocess 
import tempfile

def makePatserAnnotation(sequence=None,patser_file=None,matrix=None,seqname="NoName_sequence",scorecut=None):
    """
    Class method for generating annotation tracks of patser hits
    """
    annot = None
    
    if (sequence and not patser_file):
        annot = _patserAnnotSeq(sequence, matrix, seqname,scorecut=scorecut)
    elif (not sequence and patser_file):
        annot = _parsePatser(patser_file,matrix,scorecut=scorecut)
    else:
        annot = _parsePatser(patser_file,matrix,scorecut=scorecut)
        annot.setParentSeq(sequence)
    
    return annot

def _patserAnnotSeq(sequence,matrixfile,name,scorecut=None):

    # make sequence patser-readable *cough* BULLSHIT *cough*.  Honestly I don't really know why this offends me so much.  Oh for a world where everything just read fasta files...
    seq = name + " \\" + sequence + "\\"

    # build arg list for patser
    patserArgs = [os.environ["HOME"] + "/bin/patser",
                  "-m" + matrixfile,
                  #"-w", # weight matrix
                  #"-v", # vertical matrix
                  "-A a:t 1 c:g 1", # alphabet - no weighting for now, could calculate from sequence or from genome
                  "-c" # score complements
                 ]
    # create a tempfile for output
    stdout_temp = tempfile.NamedTemporaryFile()

    # start a subprocess for patser 
    patProc = subprocess.Popen(patserArgs,stdin=subprocess.PIPE,stdout=stdout_temp)
    
    patProc.communicate(input=seq)
    
    stdout_temp.flush()
    os.fsync(stdout_temp)
    stdout_temp.seek(0)
    
    return _parsePatser(stdout_temp,matrixfile,sequence=sequence,scorecut=scorecut)

def _parsePatser(input_fh,matrix,sequence=None,scorecut=None):
    """
    parse Patser output (contained in filehandle input_fh)
    """

    patAnnot = annotation_track.AnnotationTrack(length=len(sequence),parentSeq=sequence)
    mat_width = 0
    sugg_scorecut = 0.0
    sugg_pvalcut = 0.0
    strand = None
     
    matrix_name = re.search("\S*/(\S+?)\.",matrix).group(1)

    # regexes for extracting width and score cutoff info from file
    width_re = re.compile('\s*width of the alignment matrix:\s+(\d+)')
    scorecut_re = re.compile('\s*numerically calculated cutoff score:\s+(\S+)')
    pvalcut_re = re.compile('\s*ln\(numerically calculated cutoff p-value\):\s+(\S+)')

    # find width, score cutoff, pvalue cutoff from header
    for l in input_fh:
        #print l
        if (mat_width == 0) and width_re.match(l):
            mat_width = int(width_re.match(l).group(1))
            continue
        elif (sugg_scorecut == 0) and scorecut_re.match(l):
            sugg_scorecut = float(scorecut_re.match(l).group(1))
            continue
        elif (sugg_pvalcut == 0) and pvalcut_re.match(l):
            sugg_pvalcut = float(pvalcut_re.match(l).group(1))
            break

    # default is to set the score cutoff to the suggested value from the file
    if (not scorecut) and (scorecut != 0):
        scorecut = sugg_scorecut

    # read lines, throw out everything below the cutoff.  Create new features in the annotation track for sites passing cutoff
    for l in input_fh:
        #print l
        match = re.match("\s*(\S+)\s+position=\s+(\d+)(\w*)\s+score=\s+(\S+)(\s+ln\(p-value\)=\s+(\S+))?",l)
        
        if match and (float(match.group(4)) >= scorecut):
            #print "MATCH"
            bs_start = 0;
            bs_end = 0;
            if match.group(3) == "C": # e.g. if hit is on minus strand
                bs_start = int(match.group(2)) 
                bs_end = int(match.group(2)) + mat_width
                strand = "-"
            else:
                bs_start = int(match.group(2))
                bs_end = int(match.group(2)) + mat_width
                strand = "+"

            site = annotation_track.BindingSiteFeature(start=bs_start,end=bs_end,
                                                       tags={'matwidth' : mat_width,
                                                             'scorecut' : scorecut,
                                                             'pvalcut' : sugg_pvalcut,
                                                             'score' : float(match.group(4)),
                                                             'strand' : strand,
                                                             'motif_name' : matrix_name})
            if match.group(6):
                site.addTag(('pval',float(match.group(6))))
            patAnnot.addFeature(site)

    input_fh.close()
    return patAnnot
