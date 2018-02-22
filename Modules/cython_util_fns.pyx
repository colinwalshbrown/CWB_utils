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
uint32_DTYPE = np.uint32
ctypedef np.float_t fl_DTYPE_t
ctypedef np.uint32_t uint32_DTYPE_t
ctypedef np.str str_DTYPE_t
#char_DTYPE = np.char
#ctypedef np.char_t


# for each position x in m, add m[x] to the previous ext positions in e
def minus_extend(np.ndarray[fl_DTYPE_t] m,int alength, int ext):
    
    cdef int x = 0
    cdef int lower = 0
    cdef int ext_t = ext
    cdef long unsigned int alen = alength 
    cdef int y = 0
    cdef np.ndarray[fl_DTYPE_t] e = np.zeros(alength,dtype=fl_DTYPE)

    while x < alen:
        lower = x - ext
        if lower < 0:
            lower = 0
        y = lower

        if m[x] == 0:
            x += 1
            continue

        while y < x:
            e[y] += m[x]
            y += 1
#        print e[lower:x]
        x += 1

    return e

# for each position x in p, add p[x] to the next ext positions in e
def plus_extend(np.ndarray[fl_DTYPE_t] p,int alength, int ext):
    
    cdef int x = 0
    cdef int upper = 0
    cdef int ext_t = ext
    cdef long unsigned int alen = alength
    cdef int y = 0
    cdef np.ndarray[fl_DTYPE_t] e = np.zeros(alength,dtype=fl_DTYPE)

    while x < alen:
        upper = x + ext_t + 1
        if upper > alen:
            upper = alen
        y = x + 1

        if p[x] == 0:
            x += 1
            continue

        while y < upper:
            e[y] += p[x]
            y += 1
        x += 1

    return e

def patser_array_overlap(np.ndarray[str_DTYPE_t] chr_arr, np.ndarray[uint32_DTYPE_t, ndim=2] start_arr, np.ndarray[uint32_DTYPE_t, ndim=2] end_arr, char *chrom, long start, long end):

    cdef long ol_st = -1
    cdef long ol_end = -1
    cdef int idx = 0
    cdef unsigned int length = len(chr_arr)
    cdef unsigned long a_start = start_arr[0][0][0]
    cdef unsigned long a_end = end_arr[0][0][0]
    cdef char *a_chrom = chr_arr[0][0][0]
    cdef np.ndarray[np.int8_t] hit_arr = np.zeros(len(chr_arr),dtype=np.int8)

    print "hi"
    
    while (idx < length):

        a_start = start_arr[idx][0][0]
        a_end = end_arr[idx][0][0]
        a_chrom = chr_arr[idx][0][0]
        
        if (a_chrom == chrom):
            if ((a_start <= end) and (a_start >= start)):
                #ol_st = a_start
                hit_arr[idx] = 1
            elif ((start >= a_start) and (start <= a_end)):
                #ol_st = start
                hit_arr[idx] = 1
        
            if ((a_end <= end) and (a_end >= start)):
                #ol_end = a_end
                hit_arr[idx] = 1
            elif ((end <= a_end) and (end > a_start)):
                #ol_end = end
                hit_arr[idx] = 1

        idx += 1

    return hit_arr
        
    #if (ol_st > 0) and (ol_end > 0):
    #    return (ol_st,ol_end)
    #else:
    #    return (0,0)

