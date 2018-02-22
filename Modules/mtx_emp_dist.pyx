#!/usr/bin/env python
#
# cython_util_fns.py: some utility functions sped up using cython
#

import sys
import numpy as np

# compile-time numpy import
cimport numpy as np

# types for arrays
fl_DTYPE = np.float
int_DTYPE = np.int
ctypedef np.float_t fl_DTYPE_t
ctypedef np.int_t int_DTYPE_t

def calc_emp_dist(np.ndarray[fl_DTYPE_t] meta,np.ndarray[fl_DTYPE_t] mt1,np.ndarray[fl_DTYPE_t] mt2,np.ndarray[int_DTYPE_t] s1,np.ndarray[int_DTYPE_t] s2,int reps):

    cdef int i
    cdef np.ndarray[int_DTYPE_t] combseqs = np.zeros((len(s1) + len(s2)),dtype=int_DTYPE)
    cdef np.ndarray[fl_DTYPE_t] llrs = np.zeros(reps,dtype=fl_DTYPE) 
    
    combseqs[0:(len(s1) - 1)] = s1
    combseqs[len(s1):] = s2
    
    for i in range(reps):
        llrs[i] = (-2 * np.log(mtx_likelihood(meta,combseqs))) + (2 * np.log((mtx_likelihood(mt2,s2) + mtx_likelihood(mt1,s1))))

    return llrs

def mtx_likelihood(np.ndarray[fl_DTYPE_t] mtx, np.ndarray[int_DTYPE_t] seqs):

    cdef int i    
    cdef float likelihood

    for i in range(len(seqs)):
        likelihood += mtx[seqs[i]]

    return likelihood
