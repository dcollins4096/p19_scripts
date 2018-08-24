if 1:
    import matplotlib
    matplotlib.use('Agg')
    import math
    import matplotlib.pyplot as plt
    import numpy as np
    nar = np.array
    import time
    from importlib import reload
    import pyximport; pyximport.install()
    import particle_ops
    import particle_grid_mask
#from p19b_select_particle_callback import *
#from p19_tools import *
    import h5py
    import copy 
    import pdb
    from importlib import reload
    from collections import defaultdict
    import weakref
if 0:
    import loop_tools
    reload(loop_tools)
"""
trackage.
the tool package to get history for individual particles
"""


class track_manager():
    """Manages tracks.
    track_manager[field] returns a list of [particle, time]
    so 
    for i in track_manager[field]:
        plt.plot(track_manager.time,i)
    plots the field vs. time.
    self.ingest(snapshot) populates the master list
    """
    def __init__(self,my_loop):
        self.my_loop = weakref.proxy(my_loop)
        self.particle_ids=nar([],dtype='int64')
        self.core_ids = nar([],dtype='int64')
        self.frames = nar([],dtype='int64')
        self.times = nar([],dtype='float64')
        self.track_dict={}
        self.shape=(0,0) #particles, frames

    def ingest(self,snapshot):
        particle_ids = copy.copy(snapshot.target_indices)
        if snapshot.core_id not in self.core_ids:
            #this might not be the best place for the parent step.
            core_ids = np.ones_like(particle_ids) * snapshot.core_id
            self.core_ids = np.append(self.core_ids, core_ids)
            self.particle_ids = np.append(self.particle_ids, particle_ids)
        particle_start = np.where(self.particle_ids==particle_ids[0])[0][0]
        particle_end=particle_start+particle_ids.size

        if snapshot.frame not in self.frames:
            self.frames=np.append(self.frames,snapshot.frame)
            self.times=np.append(self.times,snapshot.ds['InitialTime'])

        frame_id = np.where(self.frames == snapshot.frame)[0][0]

        for field in snapshot.field_values:
            current_shape = self[field].shape
            new_shape = [self.particle_ids.size,
                         self.frames.size]
            temp_frame = np.zeros(new_shape,dtype='float64')
            old_slice = (slice(None,current_shape[0]),
                         slice(None,current_shape[1]))
            temp_frame[old_slice]= self[field]
            new_slice = (slice(particle_start,particle_end),
                         slice(frame_id,frame_id+1))
            nuggle=np.array(snapshot.field_values[field])
            nuggle.shape=(particle_ids.size,1)
            temp_frame[new_slice]=nuggle
            self[field]=temp_frame


    def p(self,particle_list,field):
        output = None
        for particle in particle_list:
            loc = np.where(self.particle_ids == particle)
            parts = self[field][loc,:]
            if output is None:
                output = parts
            else:
                output= np.append(output, parts,axis=0)
        return output
    def c(self,core_list,field):
        if type(core_list) is int:
            core_list = [core_list]
        nf = self.frames.size
        output = None
        for core in core_list:
            loc = self.core_ids == core
            core_values = self[field][loc,:]
            if output is None:
                output = core_values
            else:
                output = np.append(output,core_values,axis=0)
        return output
    def __getitem__(self,item):
        output = None
        if item in ['particles', 'particle_ids']:
            return self.particle_ids
        if item in ['core_ids','cores']:
            return self.core_ids
        if item not in self.track_dict:
            self.track_dict[item]=np.zeros(self.shape,dtype='float64')
        return self.track_dict[item]
    def __setitem__(self,item,value):
        self.track_dict[item]=value
    def setup_fields(self,field_list=None):
        for frame in self.my_loop.field_list:
            for core_id in self.my_loop.core_list:
                this_snapshot = looper.make_snapshot(frame,core_id)
                if field_list is None:
                    local_field_list = this_snapshot.field_values.keys()
                else:
                    local_field_list = field_list
                for field in local_field_list:
                    self[field].ingest(this_snapshot)


        
