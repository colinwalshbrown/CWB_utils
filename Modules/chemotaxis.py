#!/usr/bin/env python
#
#	chemotaxis.py:
#            Classes for running and analyzing chemotaxis simulations
#
#

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as anim
import time


# module-specific imports
import chemotaxis_sim_cython_fns as chem_c
AGENT_DTYPE = np.dtype([('x','f8'),
                        ('y','f8'),
                        ('ch_n_1','f8'),
                        ('ch_n_2','f8'),
                        ('ch_n_3','f8'),
                        ('heading','f8'),
                        ('atype','i8'),
                        ('tprob_center','f8'),
                        ('tprob_slope','f8'),
                        ('tangle_center','f8'),
                        ('tangle_slope','f8'),
                        ('mem_param','f8'),
                        ('speed','f8'),
                        ('sec_ampl','f8'),
                        ('sec_xspr','f8'),
                        ('sec_yspr','f8')])

class Agent(object):
   def __init__(self,atype=None,tprob_center=None,tprob_slope=None,tangle_center=None,
                 tangle_slope=None,mem_param=None,speed=None,sec_ampl=None,sec_xspr=None,sec_yspr=None):
      """
      A single agent in the simulation
      
      atype: Agent type (int)
      tprob_center: center parameter for the agent's tumble probability function
      tprob_slope: slope ('steepness' / Hill coeff) for agent's tumble probability function
      tangle_center: center parameter for tumble angle limit function
      tangle_slope: slope ('steepness' / Hill coeff) for agent's tumble probability function
      speed: distance agent moves in a single run
      sec_ampl: amplitude of the agent's chemoattractant secretion function
      sec_xspr: x-spread for the agent's chemoattractant secretion function
      sec_yspr: y-spread for the agent's chemoattractant secretion function
      """

      self.atype = atype
      self.tprob_center = tprob_center
      self.tprob_slope = tprob_slope
      self.tangle_center = tangle_center
      self.tangle_slope = tangle_slope
      self.mem_param = mem_param
      self.speed = speed
      self.sec_ampl = sec_ampl
      self.sec_xspr = sec_xspr
      self.sec_yspr = sec_yspr
      self.x = 0.0
      self.y = 0.0
      self.ch_n_1 = 0.0
      self.ch_n_2 = 0.0
      self.ch_n_3 = 0.0
      self.heading = 0.0

   def __str__(self):
      p = """
          ********************
          atype = %d\t\t\t#    x = %f
          tprob_center = %f\t#    y = %f
          tprob_slope = %f\t#    ch_n_1 = %f
          tangle_center = %f\t#    ch_n_2 = %f
          tangle_slope = %f\t#    ch_n_3 = %f
          mem_param = %f\t#    heading = %f 
          speed = %f\t\t#    
          sec_ampl = %f\t\t#
          sec_xspr = %f\t\t#
          sec_yspr = %f\t\t#
          """ % (self.atype, self.x, self.tprob_center, self.y, self.tprob_slope, self.ch_n_1, self.tangle_center, self.ch_n_2,
                 self.tangle_slope, self.ch_n_3, self.mem_param, self.heading, self.speed, self.sec_ampl, self.sec_xspr, self.sec_yspr)
      return p

   def get_dtype(self):
      self_arr = np.ndarray(1,dtype=AGENT_DTYPE)
      self_arr[0]['x'] = self.x
      self_arr[0]['y'] = self.y
      self_arr[0]['heading'] = self.heading
      self_arr[0]['ch_n_1'] = self.ch_n_1
      self_arr[0]['ch_n_2'] = self.ch_n_2
      self_arr[0]['ch_n_3'] = self.ch_n_3
      self_arr[0]['atype'] = self.atype
      self_arr[0]['tprob_center'] = self.tprob_center
      self_arr[0]['tprob_slope'] = self.tprob_slope
      self_arr[0]['tangle_center'] = self.tangle_center
      self_arr[0]['tangle_slope'] = self.tangle_slope
      self_arr[0]['mem_param'] = self.mem_param
      self_arr[0]['speed'] = self.speed
      self_arr[0]['sec_ampl'] = self.sec_ampl
      self_arr[0]['sec_xspr'] = self.sec_xspr
      self_arr[0]['sec_yspr'] = self.sec_yspr

      return self_arr

          
          
          
          
"""
   def init_start_loc(self,arena_size):
      self.x = np.random.uniform(low=0.00000001,high=arena_size[0])
      self.y = np.random.uniform(low=0.00000001,high=arena_size[1])
      self.heading = np.random.uniform(low=0.0,high=2*np.pi)
      self.ch_n_1 = arena_fn((self.x,self.y))
      self.ch_n_2 = arena_fn((self.x - arena_size[0]*0.0001,self.y - arena_size[1]*0.0001))
      self.ch_n_3 = arena_fn((self.x - arena_size[0]*0.0002,self.y - arena_size[1]*0.0002))
"""

class AgentFactory(object):

   def __init__(self,atype=None,tprob_center=None,tprob_slope=None,tangle_center=None,
                 tangle_slope=None,mem_param=None,speed=None,sec_ampl=None,sec_xspr=None,sec_yspr=None):
      """
      Factory Class for Agents

      atype: Agent type (int)
      tprob_center: center parameter for the agent's tumble probability function
      tprob_slope: slope ('steepness' / Hill coeff) for agent's tumble probability function
      tangle_center: center parameter for tumble angle limit function
      tangle_slope: slope ('steepness' / Hill coeff) for agent's tumble probability function
      speed: distance agent moves in a single run
      sec_ampl: amplitude of the agent's chemoattractant secretion function
      sec_xspr: x-spread for the agent's chemoattractant secretion function
      sec_yspr: y-spread for the agent's chemoattractant secretion function
      """
      self.atype = atype
      self.tprob_center = tprob_center
      self.tprob_slope = tprob_slope
      self.tangle_center = tangle_center
      self.tangle_slope = tangle_slope
      self.mem_param = mem_param
      self.speed = speed
      self.sec_ampl = sec_ampl
      self.sec_xspr = sec_xspr
      self.sec_yspr = sec_yspr
      self.constructor = self._agent_iter()

   def get_agent(self):
      """
      return a single agent
      """
      return Agent(self.atype,self.tprob_center,self.tprob_slope,self.tangle_center,
                   self.tangle_slope,self.mem_param,self.speed,self.sec_ampl,self.sec_xspr,self.sec_yspr)

   def _agent_iter(self):
      """
      iterator object for agents of the given type - might add some checks later?
      """
      while 1:
         yield Agent(self.atype,self.tprob_center,self.tprob_slope,self.tangle_center,
                     self.tangle_slope,self.mem_param,self.speed,self.sec_ampl,self.sec_xspr,self.sec_yspr)

   def next(self):
      return self.constructor.next()

                 
class AgentCollection(object):

   def __init__(self,st_agents=None):
      """
      AgentCollection(st_agents=None):
      
      a (mutable) agent collection - assemble agents here, then
      call get_agent_array() to return a

      st_agents: list of agents to initialize the collection
      """
      if st_agents:
         self.agents = st_agents
         self.length = len(st_agents)
      else:
         self.agents = []
         self.length = 0

   def add_agents(self,n,agent_factory):
      """
      add n agents to the collection using AgentFactory object agent_factory
      """
      for i in range(n):
         self.agents.append(agent_factory.next())

   def add_agent_obj(self,agent):
      """
      Add a single Agent object to the collection
      """
      self.agents.append(agent)

   def init_all(self,arena_size,arena_fn):
      """
      set random start positions and sensor values for all agents in collection
      """
      for agt in self.agents:
               agt.x = np.random.uniform(low=0.00000001,high=arena_size[0])
               agt.y = np.random.uniform(low=0.00000001,high=arena_size[1])
               agt.heading = np.random.uniform(low=0.0,high=2*np.pi)
      agt_arr = self.get_agent_array()
      arena_fn = chem_c.get_gauss_fn(agt_arr)
      for agt in self.agents:
         agt.ch_n_1 = arena_fn.evaluate(agt.x,agt.y)
         agt.ch_n_2 = arena_fn.evaluate(agt.x - arena_size[0]*0.0001,agt.y - arena_size[1]*0.0001)
         agt.ch_n_3 = arena_fn.evaluate(agt.x - arena_size[0]*0.0002,agt.y - arena_size[1]*0.0002)

   def build_agent_coord_array(self,coord_array):
      ar = np.ndarray(len(self.agents),dtype=AGENT_DTYPE)
      for (i,a) in enumerate(self.agents):
         ar[i]['x'] = coord_array[i]['x']
         ar[i]['y'] = coord_array[i]['y']
         ar[i]['heading'] = coord_array[i]['heading']
         ar[i]['ch_n_1'] = coord_array[i]['ch_n_1']
         ar[i]['ch_n_2'] = coord_array[i]['ch_n_2']
         ar[i]['ch_n_3'] = coord_array[i]['ch_n_3']
         ar[i]['atype'] = a.atype
         ar[i]['tprob_center'] = a.tprob_center
         ar[i]['tprob_slope'] = a.tprob_slope
         ar[i]['tangle_center'] = a.tangle_center
         ar[i]['tangle_slope'] = a.tangle_slope
         ar[i]['mem_param'] = a.mem_param
         ar[i]['speed'] = a.speed
         ar[i]['sec_ampl'] = a.sec_ampl
         ar[i]['sec_xspr'] = a.sec_xspr
         ar[i]['sec_yspr'] = a.sec_yspr
      return ar

   def get_agent_array(self):
      """
      return a (static) numpy ndarray of AGENT_DTYPE agents from the (mutable) AgentCollection object
      """
      ar = np.ndarray(len(self.agents),dtype=AGENT_DTYPE)
      for (i,a) in enumerate(self.agents):
         ar[i]['x'] = a.x
         ar[i]['y'] = a.y
         ar[i]['heading'] = a.heading
         ar[i]['ch_n_1'] = a.ch_n_1
         ar[i]['ch_n_2'] = a.ch_n_2
         ar[i]['ch_n_3'] = a.ch_n_3
         ar[i]['atype'] = a.atype
         ar[i]['tprob_center'] = a.tprob_center
         ar[i]['tprob_slope'] = a.tprob_slope
         ar[i]['tangle_center'] = a.tangle_center
         ar[i]['tangle_slope'] = a.tangle_slope
         ar[i]['mem_param'] = a.mem_param
         ar[i]['speed'] = a.speed
         ar[i]['sec_ampl'] = a.sec_ampl
         ar[i]['sec_xspr'] = a.sec_xspr
         ar[i]['sec_yspr'] = a.sec_yspr
      return ar

   def list_agents_by_type(self):
      """
      return a list of arrays, where each array is the agent of a single type (ID int)
      and index == agent ID
      """
      type_list = []
      for agt in self.agents:
         if agt.atype <= len(type_list):
            type_list[agt.atype - 1].append(agt.get_dtype())
         else:
            type_list.append([agt.get_dtype(),])

      return [np.array(a,dtype=AGENT_DTYPE) for a in type_list]

   def __len__(self):
      return len(self.agents)

class ChemotaxisSimulation(object):

   def __init__(self,arena_fn_gen,agents,arena_size=None):
      # initialize variables from args
      
      self.arena_fn_gen = arena_fn_gen # generator function for arena_fn functions at each timestep
      if not agents:
         self.agents = AgentCollection() # CHANGE AFTER WRITING AGENTCOLLECTION
      else:
         self.agents = agents # AgentCollection object
        
      if not arena_size:
         self.arena_size = (500,500)
      else:
         self.arena_size = arena_size

      self.nagents = len(self.agents)
      self.agents.init_all(self.arena_size,self.arena_fn_gen(self.agents.get_agent_array()))

      self.cur_sim = None
      self.cur_anim = None

   def run_chemotaxis_sim(self,nsteps=1000,reinitialize_agents=True,animation=None):
      """
      run the simulation using the specified set of agents

      parameters:
          nsteps = number of simulation steps, def. 1000
          reinitialize_agents = reset agent positions before running simulation, def. True
          animation = file name to save animation; if None no file saved
      """
      # fill agent array
      if reinitialize_agents:
         self.agents.init_all(self.arena_size,self.arena_fn_gen(self.agents.get_agent_array()))
      agent_array = self.agents.get_agent_array()
    
      sim_mtx = chem_c.build_sim_mtx(nsteps,self.nagents,agent_array,self.arena_size[0],self.arena_size[1])
      self.cur_anim = None
      if animation:
         self.animate(animation)
    
      self.cur_sim = sim_mtx

   def snapshot(self,timestep,outfile,format="jpeg"):
      if not self.cur_anim:
         self.cur_anim = AnimatedChemotaxisSim(self)
      snap_fig = self.cur_anim.snapshot(timestep)
      if outfile:
         snap_fig.savefig(outfile,format=format)

      return snap_fig
      
   def animate(self,anim_file,anim_interval=20,anim_blit=True,px_size=1):
      if not self.cur_anim:
         self.cur_anim = AnimatedChemotaxisSim(self)
      self.cur_anim.px_size=px_size
      self.cur_anim.render_animation(anim_file,anim_interval,anim_blit)
      #animation.save(anim_file, fps=30, extra_args=['-vcodec', 'libx264'])

class AnimatedChemotaxisSim(object):
   """
   Compute & render an animated visualization of a chemotaxis simulation
   """
    
   def __init__(self,
                sim,
                arena_vmax=None,
                fig_size=None,
                arena_cmap='Greys',
                markers=['.','o','v','^','<','>','s','+','x','D','d'],
                colors=['r','g','b','c','m','y','k'],
                line_styles=['-','--','-.',':'],
                px_size = 1,
                show_tracks=False,
                show_ends=True,
                show_starts=False):
      """
      sim: ChemotaxisSimulation object
      arena_vmax: ceiling value passed to plt.imshow()
      arena_pal: matplotlib color pallet to use for arena plotting
      fig_size: output size for the animation
      markers: markers for each agent type
      colors: colors for agents & traces
      line_styles: line styles for agent traces
      """
      self.sim = sim
      self.agent_types = len(sim.agents.list_agents_by_type())
      self.agent_array = sim.agents.get_agent_array()
      if sim.cur_sim == None:
         sim.run_chemotaxis_sim()
      self.sim_mtx = sim.cur_sim
      self.arena_size = sim.arena_size
      self.px_size = px_size
                    
      if not arena_vmax:
         self.arena_vmax = np.max(chem_c.build_arena(chem_c.get_gauss_fn(self.agent_array),self.arena_size[0],self.arena_size[1],px_size)) * 2
      else:
         self.arena_vmax = arena_vmax
            
      if not fig_size:
         self.figsize = (15,15)
      else:
         self.figsize = fig_size

      self.markers = markers
      self.colors = colors
      self.line_styles = line_styles

      #self.build_arena_cfn = build_arena_cfn
            
            
      # initialize matplotlib figure and axis
      self.fig = plt.figure(figsize=self.figsize)
      self.ax = self.fig.add_subplot(111, aspect='equal', autoscale_on=False,
                                     xlim=(0,self.arena_size[0]), ylim=(0,self.arena_size[1]))
      # initialize style pickers:
      self.start_markers = self.get_agent_type_plot_style(self.agent_types,colors,markers)
      self.end_markers = self.get_agent_type_plot_style(self.agent_types,colors,markers[:self.agent_types])
      self.lines = self.get_agent_type_plot_style(self.agent_types,colors,line_styles)
      """
      for i in range(self.agent_types):
         try:
            a_col = self.colors_iter.next()
         except StopIteration:
            self.colors_iter = self.get_plot_color()
            a_col = self.colors_iter.next()

         try:
            a_stmk = self.marker_iter.next()
         except StopIteration:
            self.marker_iter = self.get_plot_marker()
            a_stmk = self.marker_iter.next()

         try:
            a_endmk = self.marker_iter.next()
         except StopIteration:
            self.marker_iter = self.get_plot_marker()
            a_endmk = self.marker_iter.next()

         try:
            a_lin = self.linestyles_iter.next()
         except StopIteration:
            self.linestyles_iter = self.get_line_style()
            a_lin = self.linestyles_iter.next()
            
         self.start_markers.append(a_col + a_stmk)
         self.end_markers.append(a_col + a_endmk)
         self.lines.append(a_col + a_lin)
      """
      # initialize drawables:
      #
      # arena - background plot of attractant concentration (static or can be re-calculated each frame)
      self.empty_arena = np.zeros(self.arena_size)
      self.arena = self.ax.imshow(self.empty_arena,cmap=arena_cmap,vmax=self.arena_vmax) # vmax will default to 0 if imshow initialized with zeros
        
      # starting location of the agents - static for all frames
      self.agent_starts = []
      if show_starts:
         self.agent_starts = [self.ax.plot(self.agent_array[i]['x'],self.agent_array[i]['y'],self.start_markers[self.agent_array[i]['atype'] - 1],markersize=5)
                                       for i in range(len(self.agent_array))]
      self.agent_tracks = []
      if show_tracks:
         self.agent_tracks = [self.ax.plot([],[],self.lines[self.agent_array[x]['atype'] - 1],alpha=0.25)[0] for x in range(len(self.agent_array))]

      self.show_ends = []
      if show_ends:
         self.agent_ends = [self.ax.plot([],[],self.end_markers[self.agent_array[x]['atype'] - 1],markersize=5)[0] for x in range(len(self.agent_array))]

      """
      # Should be able to skip this w/ C implementation of build_arena
      # build a MxN matrix template, where each cell = (i,j), i.e. a tuple of the cell's own coordinates
      self.arena_template = np.zeros(self.arena_size,dtype=np.dtype([('x','f8'),('y','f8')]))
      self.arena_template[:,:]['x'] = np.vstack([np.linspace(1,self.arena_size[0],num=self.arena_size[0])] * self.arena_size[1])
      self.arena_template[:,:]['y'] = np.vstack([np.linspace(1,self.arena_size[1],num=self.arena_size[1])] * self.arena_size[0]).T
      """

   def style_iter(self,style_list):
      for s in style_list:
         yield s

   def get_agent_type_plot_style(self,n_agent_types,colors,styles):
      colors_iter = self.style_iter(colors)
      styles_iter = self.style_iter(styles)
      plot_styles = []
      
      for i in range(n_agent_types):
         try:
            a_col = colors_iter.next()
         except StopIteration:
            colors_iter = self.style_iter(colors)
            a_col = colors_iter.next()

         try:
            a_style = styles_iter.next()
         except StopIteration:
            self.marker_iter = self.style_iter(styles)
            a_style = self.styles_iter.next()

         plot_styles.append(a_col + a_style)

      return plot_styles

   def _anim_init(self):
      """
      init method for anim.FuncAnimation
      """
        
      # clear the arena
      self.arena.set_array(self.empty_arena)
       
      # clear agent tracks and endpoints
      for i in range(self.agent_array.shape[0]):
         if len(self.agent_tracks) > 0:
            self.agent_tracks[i].set_data([],[])
         if len(self.agent_ends) > 0:
            self.agent_ends[i].set_data([],[])
            
      all_plots = self.agent_tracks + self.agent_ends + [self.arena]
      return all_plots # remember for the future that returned drawables need to be in the form of an iterable!

   def _animate(self,i,diag=True,diag_check=10):
      """
      animate method for anim.FuncAnimation
      """
      diagnostics = False
      if diag and i % diag_check == 0:
         diagnostics = True
      # start a timer
      if diagnostics:
         print "animating step %d..." % (i,)
      st_time = time.time()

      # build agent array for this timestep
      t_array = self.sim.agents.build_agent_coord_array(self.sim_mtx[:,i])
      t_array_time = time.time() - st_time
        
      # build the arena for this timestep
      arena_fn = chem_c.get_gauss_fn(t_array)
      arena_fn_time = time.time() - (t_array_time + st_time)
      self.arena.set_array(chem_c.build_arena(arena_fn,self.arena_size[0],self.arena_size[1],self.px_size).T)
        
      # arena build time
      arena_build_time = time.time() - (arena_fn_time + t_array_time + st_time)
        
      if diagnostics:
         #print "\t MAX: %f, MIN: %f" % (np.max(a),np.min(a))
         
         print "\tt_array build time: %f" % (t_array_time,)
         print "\tarena_fn build time: %f" % (arena_fn_time,)
         print "\tarena build time: %f" % (arena_build_time,)
        
      # plot track and endpoint for each agent at timestep i
      for j in range(self.sim_mtx.shape[0]):
         if len(self.agent_tracks) > 0:
            self.agent_tracks[j].set_data(self.sim_mtx[j,:i]['x'],self.sim_mtx[j,:i]['y'])
         if len(self.agent_ends) > 0:
            self.agent_ends[j].set_data(self.sim_mtx[j,i]['x'],self.sim_mtx[j,i]['y'])
        
      agent_plot_time = time.time() - (st_time + arena_build_time)
      if diagnostics:
         print "\tagent plot time: %f" % (agent_plot_time,)
        
      all_plots = self.agent_tracks + self.agent_ends + [self.arena]
        
      tot_time = time.time() - st_time
      if diagnostics:   
         print "done! total frame time: %f" % (tot_time,)
        
      return all_plots

   def render_animation(self,anim_file,anim_fps=30,anim_interval=20,anim_blit=True):
      interval_frames = self.sim_mtx[::anim_interval]
      print "Animating %d steps..." % (anim_frames)
      # render & save animation
      animation = anim.FuncAnimation(self.fig, self._animate, frames=anim_frames,
                                     blit=anim_blit, init_func=self._anim_init)
        
      animation.save(anim_file, fps=anim_fps, extra_args=['-vcodec', 'libx264'])

   def snapshot(self,
                timestep,
                fig_size=None,
                arena_cmap='Greys',
                markers=['.','o','v','^','<','>','s','+','x','D','d'],
                colors=['r','g','b','c','m','y','k'],
                line_styles=['-','--','-.',':'],
                show_arena=True,
                show_tracks=False,
                show_ends=True,
                show_starts=False):

      start_markers = self.get_agent_type_plot_style(self.agent_types,colors,markers)
      end_markers = self.get_agent_type_plot_style(self.agent_types,colors,markers[self.agent_types:])
      lines = self.get_agent_type_plot_style(self.agent_types,colors,line_styles)

      print start_markers
      print end_markers
      print lines

      t_array = self.sim.agents.build_agent_coord_array(self.sim_mtx[:,timestep])
      arena_fn = chem_c.get_gauss_fn(t_array)
      arena = chem_c.build_arena(arena_fn,self.arena_size[0],self.arena_size[1],self.px_size).T
      if show_arena:
         plt.imshow(arena,interpolation='none',cmap=arena_cmap)
      for j in range(self.agent_array.shape[0]):
         if show_starts:
            plt.plot(self.sim_mtx[j,0]['x'],self.sim_mtx[j,0]['y'],lines[self.agent_array[j]['atype'] - 1])
         if show_tracks:
            plt.plot(self.sim_mtx[j,:timestep]['x'],self.sim_mtx[j,:timestep]['y'],lines[self.agent_array[j]['atype'] - 1])
         if show_ends:
            plt.plot(self.sim_mtx[j,timestep]['x'],self.sim_mtx[j,timestep]['y'],end_markers[self.agent_array[j]['atype'] - 1])
      return plt.gcf()


# SCRAAAAAPS
      
"""

     fig = plt.figure(figsize=self.figsize)
      ax = self.fig.add_subplot(111, aspect='equal', autoscale_on=False,
                                xlim=self.arena_xlim, ylim=self.arena_ylim)
      if not anim_frames:
         self.anim_frames = self.agent_array.shape[1]
        
      # render & save animation
      animation = anim.FuncAnimation(fig, self._animate, frames=anim_frames,
                                     interval=anim_interval, blit=anim_blit, init_func=self._anim_init)


   def build_arena(arena_size,arena_fn):
      arena_template = np.zeros(arena_size,dtype=np.dtype([('x','f8'),('y','f8')]))
      arena_template[:,:]['x'] = np.vstack([np.linspace(1,arena_size[0],num=arena_size[0])] * arena_size[1])
      arena_template[:,:]['y'] = np.vstack([np.linspace(1,arena_size[1],num=arena_size[1])] * arena_size[0]).T
      arena_fn_vect = np.vectorize(arena_fn)
      arena = np.zeros(arena_size)
      arena[:,:] = arena_fn_vect(arena_template)
      return arena
    
   def plot_arena(arena_size,arena_fn,subplot,agent_matrix):
      arena_template = np.zeros(arena_size,dtype=np.dtype([('x','f8'),('y','f8')]))
      arena_template[:,:]['x'] = np.vstack([np.linspace(1,arena_size[0],num=arena_size[0])] * arena_size[1])
      arena_template[:,:]['y'] = np.vstack([np.linspace(1,arena_size[1],num=arena_size[1])] * arena_size[0]).T
      arena_fn_vect = np.vectorize(arena_fn)
      arena = np.zeros(arena_size)
      arena[:,:] = arena_fn_vect(arena_template)
      plt.xlim(0,arena_size[0])
      plt.ylim(0,arena_size[1])
      plt.subplot(*subplot)
      plt.imshow(arena)
      for i in range(agent_matrix.shape[0]):
         plt.plot(agent_matrix[i,:]['x'],agent_matrix[i,:]['y'],'w',alpha=0.25)
         plt.plot(agent_matrix[i,0]['x'],agent_matrix[i,0]['y'],'r>',markersize=5)
         plt.plot(agent_matrix[i,-1]['x'],agent_matrix[i,-1]['y'],'g>',markersize=5)
"""
"""        
   def build_arena(self,arena_fn):
      
      The value of the given arena function, where arena_size = (M,N)
      
      arena_fn_vect = np.vectorize(arena_fn)
      arena = np.empty(self.arena_size)
      arena[:,:] = arena_fn_vect(self.arena_template)

      return arena
    
   def build_arena_c(self,timestep,arena_cfn):
      
      if build_arena_cfn is defined, set up agent parameter arrays and pass to the arena function
   
      x_positions = np.hstack(np.array(self.agent_array[:,timestep]['x'],dtype=np.float64))
      y_positions = np.hstack(np.array(self.agent_array[:,timestep]['y'],dtype=np.float64))
      amplitudes = np.zeros(x_positions.shape[0],dtype=np.float64) + 10
      x_spreads = np.zeros(x_positions.shape[0],dtype=np.float64) + 5
      y_spreads = np.zeros(x_positions.shape[0],dtype=np.float64) + 7.5
      arena = arena_cfn(self.arena_size[0],self.arena_size[1],x_positions,y_positions,amplitudes,x_spreads,y_spreads)
      return arena.T
"""
"""
   def _anim_init(self):
        \"""
        init method for anim.FuncAnimation
        \"""
        
        # clear the arena
        self.arena.set_array(self.empty_arena)
        
        # clear agent tracks and endpoints
        for i in range(self.agent_array.shape[0]):
            self.agent_tracks[i].set_data([],[])
            self.agent_ends[i].set_data([],[])
            
        all_plots = self.agent_tracks + self.agent_ends + [self.arena]
        return all_plots # remember for the future that returned drawables need to be in the form of an iterable!

   def _animate(self,i,diag=True,diag_check=10):
      \"""
      animate method for anim.FuncAnimation
      \"""
      diagnostics = False
      if diag and i % diag_check == 0:
         diagnostics = True
      # start a timer
      if diagnostics:
         print "animating step %d..." % (i,)
      st_time = time.time()
        
      # build the arena for this timestep
      a = None
      if self.build_arena_cfn:
         a = self.build_arena_c(i,self.build_arena_cfn)
      else:
         fn = self.fn_array[i]
         a = self.build_arena(fn)
      self.arena.set_array(a)
        
      # arena build time
      arena_build_time = time.time() - st_time
        
      if diagnostics:
         print "\t MAX: %f, MIN: %f" % (np.max(a),np.min(a))
         print "\tarena build time: %f" % (arena_build_time,),
        
      # plot track and endpoint for each agent at timestep i
      for j in range(self.agent_array.shape[0]):
         self.agent_tracks[j].set_data(self.agent_array[j,:i]['x'],self.agent_array[j,:i]['y'])
         self.agent_ends[j].set_data(self.agent_array[j,i]['x'],self.agent_array[j,i]['y'])
        
      agent_plot_time = time.time() - (st_time + arena_build_time)
      if diagnostics:
         print "\tagent plot time: %f" % (agent_plot_time,)
        
      all_plots = self.agent_tracks + self.agent_ends + [self.arena]
        
      tot_time = time.time() - st_time
      if diagnostics:   
         print "done! total frame time: %f" % (tot_time,)
        
      return all_plots
"""
