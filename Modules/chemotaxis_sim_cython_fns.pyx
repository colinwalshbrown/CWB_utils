#!/usr/bin/env python
#
# chemotaxis_sim_cython_fns.pyx: some utility functions sped up using cython for chemotaxis_sim scripts
#

import sys

import numpy as np
from time import clock

# compile-time numpy import
cimport numpy as np

# types for arrays
fl_DTYPE = np.float
uint32_DTYPE = np.uint32
ctypedef np.float_t fl_DTYPE_t
ctypedef np.uint32_t uint32_DTYPE_t
ctypedef np.float64_t dbl_DTYPE_t
#ctypedef np.str str_DTYPE_t
#char_DTYPE = np.char
#ctypedef np.char_t

# try using GSL RNG instead of numpy
"""
cdef extern from "time.h":
    ctypedef struct clock_t:
        pass
    clock_t clock()
"""
cdef extern from "gsl/gsl_rng.h":
    ctypedef struct gsl_rng_type:
        pass
    ctypedef struct gsl_rng:
        pass
    gsl_rng_type *gsl_rng_mt19937
    gsl_rng *gsl_rng_alloc(gsl_rng_type * T)
    double gsl_rng_uniform(gsl_rng * r)
    double gsl_rng_uniform_pos(gsl_rng * r)
 
cdef extern from "gsl/gsl_randist.h":
    double gamma "gsl_ran_gamma"(gsl_rng * r,double,double)
    double gaussian "gsl_ran_gaussian"(gsl_rng * r,double)
 
cdef gsl_rng *r = gsl_rng_alloc(gsl_rng_mt19937)

# external & convenience functions
cdef extern from "math.h":
    double log(double x)
    double exp(double x)
    double sin(double x)
    double cos(double x)
    double M_PI

cdef double square(double x):
    return x*x

cdef double average(long n, double[:] ls):
    cdef double s
    for k in ls:
        s += k
    return k/n

# structs for agent information
#
# coordinates:
cdef struct agent_coords:
     double x
     double y
     double ch_n_1
     double ch_n_2
     double ch_n_3
     double heading

# parameters:
cdef struct agent_params:
     unsigned int atype
     double tprob_center
     double tprob_slope
     double tangle_center
     double tangle_slope
     double mem_param
     double speed
     double sec_ampl
     double sec_xspr
     double sec_yspr

# numpy dtypes matching structs above

AGENT_COORDS_DTYPE = np.dtype({'names':['x','y','ch_n_1','ch_n_2','ch_n_3','heading'],
                               'formats':[np.float64,np.float64,np.float64,np.float64,np.float64,np.float64]},
                               align=True)

"""

    [('x',np.float64),
                               ('y',np.float64),
                               ('ch_n_1',np.float64),
                               ('ch_n_2',np.float64),
                               ('ch_n_3',np.float64),
                               ('heading',np.float64)])
"""

AGENT_PARAMS_DTYPE = np.dtype({'names':['atype','tprob_center','tprob_slope','tangle_center','tangle_slope','mem_param','speed','sec_ampl','sec_xspr','sec_yspr'],
                               'formats':[np.uint32,np.float64,np.float64,np.float64,np.float64,np.float64,np.float64,np.float64,np.float64,np.float64]},
                               align=True)
"""
    ('atype',np.uint32_t),
                               ('speed',np.float64_t),
                               ('tprob_center',np.float64_t),
                               ('tprob_slope',np.float64_t),
                               ('tangle_center',np.float64_t),
                               ('tangle_slope',np.float64_t),
                               ('sec_ampl',np.float64_t),
                               ('sec_xspr',np.float64_t),
                               ('sec_yspr',np.float64_t)])
"""
                               

# class definitions to allow callbacks for arbitrary Arena and Tumble probability / angle functions
cdef class ArenaFunction:
    #cpdef __init__(self,agent_params[:] agent_param_arr):
    """
    generic arena function base class; calculate value for a point in the arena
    """
    cdef public agent_params[:] agent_param_arr 
    cdef public agent_coords[:] agent_coord_arr
    cdef public long nagents

    def __cinit__(self,agent_params[:] agent_param_arr, agent_coords[:] agent_coord_arr, long nagents):
        self.agent_param_arr = agent_param_arr
        self.agent_coord_arr = agent_coord_arr
        self.nagents = nagents

    cpdef double evaluate(self, double x, double y)  except *:
        return 0

cdef class GaussFn(ArenaFunction):
    """
    arena function calculated as sum of gaussians across all agents; parameters: amplitude, x & y spread
    """
    cpdef double evaluate(self, double x, double y) except *:
        cdef double result = 0

        for i in range(0,self.nagents):
            try:
                result += (self.agent_param_arr[i].sec_ampl*exp(-((square((x-self.agent_coord_arr[i].x))/square(2*self.agent_param_arr[i].sec_xspr))+(square((y-self.agent_coord_arr[i].y))/square(2*self.agent_param_arr[i].sec_yspr)))))
            except Exception as E:
                print "error at agent %s, %s" % (str(self.agent_param_arr[i]),str(self.agent_coord_arr[i]))
                print str(E)
                sys.exit(0)

        return result


def get_gauss_fn(np.ndarray agent_arr):
    cdef int nagents = agent_arr.shape[0]
    cdef np.ndarray agent_coords_arr_np = np.ndarray((nagents,),dtype=AGENT_COORDS_DTYPE)
    cdef np.ndarray agent_params_arr_np = np.ndarray((nagents,),dtype=AGENT_PARAMS_DTYPE)
    cdef agent_coords [:] agent_coords_arr = agent_coords_arr_np
    cdef agent_params [:] agent_params_arr = agent_params_arr_np
    for i in range(0,nagents):
        # Initial positions & sensor values
        agent_coords_arr[i].x = agent_arr[i]['x']
        agent_coords_arr[i].y = agent_arr[i]['y']
        agent_coords_arr[i].ch_n_1 = agent_arr[i]['ch_n_1']
        agent_coords_arr[i].ch_n_2 = agent_arr[i]['ch_n_2']
        agent_coords_arr[i].ch_n_3 = agent_arr[i]['ch_n_3']
        agent_coords_arr[i].heading = agent_arr[i]['heading']
        # Single agent parameters
        agent_params_arr[i].atype = agent_arr[i]['atype']
        agent_params_arr[i].tprob_center = agent_arr[i]['tprob_center']
        agent_params_arr[i].tprob_slope = agent_arr[i]['tprob_slope']
        agent_params_arr[i].tangle_center = agent_arr[i]['tangle_center']
        agent_params_arr[i].tangle_slope = agent_arr[i]['tangle_slope']
        agent_params_arr[i].mem_param = agent_arr[i]['mem_param']
        agent_params_arr[i].speed = agent_arr[i]['speed']
        agent_params_arr[i].sec_ampl = agent_arr[i]['sec_ampl']
        agent_params_arr[i].sec_xspr = agent_arr[i]['sec_xspr']
        agent_params_arr[i].sec_yspr = agent_arr[i]['sec_yspr']
    return GaussFn(agent_params_arr,agent_coords_arr,nagents)


cdef class TumbleProbFn:
    """
    Base class for functions that relate sensor values (ratios of ch over timepoints) to tumble probabilities
    (either angle spread or tumble/run prob)
    """
    cpdef double evaluate(self, double ch, double center, double slope):
        return 0

cdef class TumbleProbLogisticFn(TumbleProbFn):
    """
    tumble/run probability modeled as a logistic function of ch
    """
    cpdef double evaluate(self, double ch, double center, double slope):
        return (0.1 + 0.8/(1 + exp(-slope*(-ch - center))))

cdef class AngleProbLogisticFn(TumbleProbFn):
    """
    angle change (turn) probability modeled as logistic function of ch
    """
    cpdef double evaluate(self, double ch, double center, double slope):
        return (M_PI/(1 + exp(slope*(-ch - center))))

cdef double tprob_logistic(double ch, double center, double slope):
    return (1/(1 + exp(-slope*(-ch - center))))

cdef double tangle_logistic(double ch, double center, double slope):
    return (M_PI/(1 + exp(slope*(-ch - center))))

cdef double mem_log_ratio_average(double ch_n, double ch_n_1, double ch_n_2, double ch_n_3, double mem_param):
    cdef double[2] ch_ratios
    ch_ratios[0] = (ch_n_1/ch_n_2)
    ch_ratios[1] = (ch_n_2/ch_n_3) # ratios for previous two timesteps
    return log((ch_n/ch_n_1)) - mem_param*log(average(2,ch_ratios)) # log-transformed ratio of current reading to average of previous 3

cdef double mem_fixed_mlt_log_avg(double ch_n, double ch_n_1, double ch_n_2, double ch_n_3, double mem_param):
    cdef double[2] ch_ratios
    ch_ratios[0] = (ch_n_1/ch_n_2) * mem_param
    ch_ratios[1] = (ch_n_2/ch_n_3) * (mem_param ** 2) # ratios for previous two timesteps
    return log((ch_n/ch_n_1) / average(2,ch_ratios)) # log-transformed ratio of current reading to average of previous 3

def build_sim_mtx(int nsteps, int nagents, np.ndarray agent_arr, double arena_x, double arena_y):
    """
    Setup and run a single simulation given an agent array, arena function class (NOT a pre-calc'd arena fn), and arena x & y bounds

    Called from ChemotaxisSimulation object
    """
    cdef np.ndarray agent_coord_time_mtx_nd = np.ndarray((nagents,nsteps),dtype=AGENT_COORDS_DTYPE)
    cdef np.ndarray agent_param_arr_nd = np.ndarray((nagents,),dtype=AGENT_PARAMS_DTYPE)
    cdef agent_coords [:,:] agent_coord_time_mtx = agent_coord_time_mtx_nd
    cdef agent_params [:] agent_param_arr = agent_param_arr_nd

    # Fill struct arrays
    for i in range(0,nagents):
        # Initial positions & sensor values
        agent_coord_time_mtx[i][0].x = agent_arr[i]['x']
        agent_coord_time_mtx[i][0].y = agent_arr[i]['y']
        agent_coord_time_mtx[i][0].ch_n_1 = agent_arr[i]['ch_n_1']
        agent_coord_time_mtx[i][0].ch_n_2 = agent_arr[i]['ch_n_2']
        agent_coord_time_mtx[i][0].ch_n_3 = agent_arr[i]['ch_n_3']
        agent_coord_time_mtx[i][0].heading = agent_arr[i]['heading']

        # Single agent parameters
        agent_param_arr[i].tprob_center = agent_arr[i]['tprob_center']
        agent_param_arr[i].tprob_slope = agent_arr[i]['tprob_slope']
        agent_param_arr[i].tangle_center = agent_arr[i]['tangle_center']
        agent_param_arr[i].tangle_slope = agent_arr[i]['tangle_slope']
        agent_param_arr[i].mem_param = agent_arr[i]['mem_param']
        agent_param_arr[i].speed = agent_arr[i]['speed']
        agent_param_arr[i].sec_ampl = agent_arr[i]['sec_ampl']
        agent_param_arr[i].sec_xspr = agent_arr[i]['sec_xspr']
        agent_param_arr[i].sec_yspr = agent_arr[i]['sec_yspr']

    # convert from c array back to a numpy ndarray for passing back to python
    agent_coord_time_mtx_nd[:,:] = run_sim(agent_coord_time_mtx, agent_param_arr, nsteps, nagents, arena_x, arena_y)
    return agent_coord_time_mtx_nd

cdef agent_coords[:,:] run_sim(agent_coords[:,:] agent_coord_mtx, agent_params[:] agent_params_arr, int nsteps, int nagents, double arena_x, double arena_y):
    cdef int i
    cdef int j
    cdef int k
    # ndarray pretty much only for initializing c array - hopefully this speeds it up?
    cdef np.ndarray step_arr_nd = np.ndarray((nagents,),dtype=AGENT_COORDS_DTYPE)
    cdef np.ndarray step_ticks_nd = np.ndarray((nagents,))
    cdef agent_coords[:] step_arr = step_arr_nd
    #cdef double[:] step_ticks = step_ticks_nd
    #cdef double st_time,end_time,ag_start,ag_end
    #tprob_fn = TumbleProbLogisticFn()
    #tangle_fn = AngleProbLogisticFn()

    for i in range(1,nsteps):
        #st_time = clock()
        for k in range(0,nagents):
            # extract array of agent_coords for timestep i
            step_arr[k] = agent_coord_mtx[k][i-1]
        # initalize ArenaFunction object for this timestep
        arena_fn = GaussFn(agent_params_arr, step_arr, nagents)
        for j in range(0,nagents):
            #ag_start = clock()
            # get coords for agent in previous timestep
            cur_coords = agent_coord_mtx[j][i-1]
            cur_params = agent_params_arr[j]
            # hopefully not too much overhead from calls for these - maybe re-code as
            # normal cdefs

            # calculate new coordinates for the agent
            agent_coord_mtx[j][i] = step_agent(arena_fn,cur_coords,cur_params,arena_x,arena_y)
            #ag_end = clock()
            #step_ticks[j] = ag_end - ag_start
        #if i % 100 == 0:
        #    end_time = clock()
        #    print "step: %f step_time: %f mean 1-agent time: %f" % (i,end_time - st_time,average(nagents,step_ticks))
    return agent_coord_mtx
    
cdef agent_coords step_agent(ArenaFunction arena_fn, agent_coords nminus1_agt, agent_params n_agt_params, double arena_x, double arena_y):
    """
    calculate a single timestep for a single agent

    this needs to be fast!  might want to swap out remaining numpy functions for c library functions
    where possible
    """
    cdef agent_coords n_agt = nminus1_agt # create a new agent (copy of previous)
    #cdef double ch_n = arena_fn.evaluate(nminus1_agt.x,nminus1_agt.y) # read chemoattractant concentration (arena_fn)
    cdef double ch_n = arena_fn.evaluate(nminus1_agt.x,nminus1_agt.y) # read chemoattractant concentration (arena_fn)
    #cdef double[2] ch_ratios
    #ch_ratios[0] = (nminus1_agt.ch_n_1/nminus1_agt.ch_n_2)
    #ch_ratios[1] = (nminus1_agt.ch_n_2/nminus1_agt.ch_n_3) # ratios for previous two timesteps
    cdef double ch_m = mem_log_ratio_average(ch_n,nminus1_agt.ch_n_1,nminus1_agt.ch_n_2,nminus1_agt.ch_n_3,n_agt_params.mem_param) # log-transformed ratio of current reading to average of previous 3
    tprob = tprob_logistic(ch_m,n_agt_params.tprob_center,n_agt_params.tprob_slope) # get a tprob using ch_m and the tprob_fn w/ agent's sensing parameters
    tprand = gsl_rng_uniform_pos(r)
    
    if tprand < tprob: # tumble (pick a new random direction using tangle_fn)
        tangle_lim = tangle_logistic(ch_m,n_agt_params.tangle_center,n_agt_params.tangle_slope)
        tarand = gaussian(r,tangle_lim) # tangle_lim sets variance of the normal random variable
        n_agt.heading += tarand
    else: # run (move one step along heading vector)
        n_agt.x += n_agt_params.speed * cos(nminus1_agt.heading) - n_agt_params.speed * sin(nminus1_agt.heading)
        n_agt.y += n_agt_params.speed * sin(nminus1_agt.heading) + n_agt_params.speed * cos(nminus1_agt.heading)
        if n_agt.x > arena_x:
            n_agt.x = 0.000000000001
        elif n_agt.x <= 0.000000000001:
            n_agt.x = arena_x
        if n_agt.y > arena_y:
            n_agt.y = 0.000000000001
        elif n_agt.y <= 0.000000000001:
            n_agt.y = arena_y

    # shift sensor values
    n_agt.ch_n_3 = n_agt.ch_n_2
    n_agt.ch_n_2 = n_agt.ch_n_1
    n_agt.ch_n_1 = ch_n
        
    return n_agt
        
cpdef double build_arena_max(ArenaFunction arena_fn,int arena_x,int arena_y):

    cdef np.ndarray arena_nd = np.zeros((arena_x,arena_y))
    cdef double[:,:] arena = arena_nd
    cdef double amax = 0
    cdef double amin = 0
    cdef int i,j
    cdef double i1,j1

    for i in range(0,arena_x):
        for j in range(0,arena_y):
            i1 = i
            j1 = j
            arena[i,j] = arena_fn.evaluate(i1+0.5,j1+0.5)
            if arena[i,j] > amax:
                amax = arena[i,j]
            elif arena[i,j] < amin:
                amin = arena[i,j]
                
    return amax


cpdef np.ndarray build_arena(ArenaFunction arena_fn,int arena_x,int arena_y,int pixelsize):

    cdef int arena_view_x = np.ceil(arena_x / pixelsize)
    cdef int arena_view_y = np.ceil(arena_y / pixelsize)
    cdef np.ndarray arena_nd = np.zeros((arena_view_x,arena_view_y),dtype=np.float64)
    cdef double[:,:] arena = arena_nd
    cdef int i,j
    cdef double i1,j1
    
    for i in range(arena_view_x):
        for j in range(arena_view_y):
            i1 = i * pixelsize
            j1 = j * pixelsize
            arena[i,j] = arena_fn.evaluate(i1+0.5,j1+0.5)

    r_arena = np.array(arena)
    return r_arena

cdef double[:,:] agent_pw_dist_array(agent_coords[:] step_agents):
    cdef np.ndarray pw_np = np.zeros((len(step_agents),len(step_agents)))
    cdef double [:] pw_dist = pw_np
    cdef double [:,:] pw_dist 

cpdef np.ndarray mean_pairwise_agt_dist(agent_coords[:,:] sim_mtx,int interval):

    cdef np.ndarray pw_np = np.zeros((sim_mtx.shape[1] / interval))
    cdef double [:] pw_dist = pw_np
    cdef agent_coords [:] sl
    cdef double [:] pr_arr
    cdef int i,j,l

    for i in np.arange(0,sim_mtx.shape[1],interval):
        sl = sim_mtx[:,i]
        l = len(sl)
        pr_arr = np.zeros(((l**2)-l,))
        for j in np.arange(l):
            pr_arr[j*l : (j*l + (l-j)) - 1] = np.sqrt((sl[j+1:].x-sl[j].x)**2 + (sl[j+1:].y - sl[j].y)**2)
        pw_dist[i] = np.mean(pr_arr)
        
    return pw_dist

cpdef np.ndarray max_grp_radius(agent_coords[:,:] sim_mtx,int interval,float radius):

    cdef np.ndarray pw_np = np.zeros((sim_mtx.shape[1] / interval))
    cdef double [:] pw_dist = pw_np
    cdef agent_coords [:] sl
    cdef 
    cdef int i,j,l
    
    for i in np.arange(0,sim_mtx.shape[1],interval):
        sl = sim_mtx[:,i]
        l = len(sl)
        agt_grps = []
        #agt_pw_dist = -np.ones(l,l))
        for j in range(l):
            agt_pw_dist = -np.ones((l,))
            for k in range(l):
                if j == k:
                    continue
                agt_pw_dist = np.sqrt((sl[j]['x']-sl[k]['x'])**2 + (sl[j]['y'] - sl[k]['y'])**2)
            agt_rad_grps = np.where((agt_pw_dist < radius) and (agt_pw_dist >=0))
            agt_grps.append(agt_rad_grps[0].shape[0])
        #rad_agts = agt_pw_dist[np.where((agt_pw_dist < radius) and (agt_pw_dist >=0))]
        #print agt_grps
        #rad_grps = [agt_pw_dist[x,:] for x in range(l)]
        #print rad_grps[0]
        #max_grp_size.append(np.max(agt_grps))
        #median_grp_size.append(np.median(agt_grps))
        #low_qrt.append(np.percentile(agt_grps,25))
        #up_qrt.append(np.percentile(agt_grps,75))
        
    return agt_grps #max_grp_size, median_grp_size, low_qrt, up_qrt


    
"""
# for each position x in m, add m[x] to the previous ext positions in e
cpdef sum_of_gaussians(np.float x, np.float y, np.ndarray[fl_DTYPE_t] x_pos, np.ndarray[fl_DTYPE_t] y_pos, np.ndarray[fl_DTYPE_t] amplitudes, np.ndarray[fl_DTYPE_t] x_sprs, np.ndarray[fl_DTYPE_t] y_sprs):

    cdef int i = 0
    cdef float v = 0
    #cdef float agent_v = 0
    cdef long unsigned int alen = x_pos.shape[0]
    #cdef int y = 0
    #cdef np.ndarray[fl_DTYPE_t] e = np.zeros(alength,dtype=fl_DTYPE)
    #print alen
    while i < alen:
        v += (amplitudes[i]*np.exp(-((np.square((x-x_pos[i]))/np.square(2*x_sprs[i]))+(np.square((y-y_pos[i]))/np.square(2*y_sprs[i])))))
        i += 1
    return v

cpdef sum_gaussians_arena(np.int arena_x, np.int arena_y, np.ndarray[fl_DTYPE_t] x_pos, np.ndarray[fl_DTYPE_t] y_pos, np.ndarray[fl_DTYPE_t] amplitudes, np.ndarray[fl_DTYPE_t] x_sprs, np.ndarray[fl_DTYPE_t] y_sprs):

    cdef unsigned int i
    cdef unsigned int j
    cdef np.ndarray[fl_DTYPE_t,ndim=2] arena = np.zeros((arena_x,arena_y),dtype=fl_DTYPE)

    for j in np.arange(0,arena_x,dtype=uint32_DTYPE):
        for i in np.arange(0,arena_y,dtype=uint32_DTYPE):
            arena[j,i] = sum_of_gaussians_c(np.float(j),np.float(i),x_pos,y_pos,amplitudes,x_sprs,y_sprs)
    return arena

cdef sum_of_gaussians_c(np.float x, np.float y, np.ndarray[fl_DTYPE_t] x_pos, np.ndarray[fl_DTYPE_t] y_pos, np.ndarray[fl_DTYPE_t] amplitudes, np.ndarray[fl_DTYPE_t] x_sprs, np.ndarray[fl_DTYPE_t] y_sprs):
    
    cdef int k = 0
    cdef float v = 0
    #cdef float agent_v = 0
    cdef long unsigned int alen = x_pos.shape[0]
    #cdef int y = 0
    #cdef np.ndarray[fl_DTYPE_t] e = np.zeros(alength,dtype=fl_DTYPE)
    #print alen
    while k < alen:
        v += (amplitudes[k]*exp(-((square((x-x_pos[k]))/square(2*x_sprs[k]))+(square((y-y_pos[k]))/square(2*y_sprs[k])))))
        k += 1
    return v
"""
