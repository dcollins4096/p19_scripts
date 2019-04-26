
import matplotlib
matplotlib.use('Agg')
import yt
import math
import matplotlib.pyplot as plt
import numpy as np
import pdb
nar = np.array
fptr = open('n_particles.txt','r')
lines=fptr.readlines()
fptr.close()
parts = np.zeros([len(lines),2])
for n,line in enumerate(lines):
    parts[n] = np.array(line.split(),dtype='int')
all_nonzero = parts[:,0][ parts[:,1] >0]
from importlib import reload

import looper
reload(looper)
import trackage
reload(trackage)
import tracks_read_write as trw
reload(trw)

if 'this_looper' not in dir():
    directory = '/scratch2/dcollins/Paper19_48/B02/u05-r4-l4-128'
    this_looper = looper.core_looper(directory= directory,
                                     sim_name = 'u05',
                                     out_prefix = 'test',
                                     target_frame = 125,
                                     frame_list =  [120], #[120,125], #[0,1,2]+list(range(10,130,10))+[125],
                                     core_list =  [96], # core_list,
                                     fields_from_grid=['x','y','z','velocity_magnitude','magnetic_field_strength',
                                                      'velocity_divergence'],
                                  )
    this_looper.get_target_indices(h5_name='u05_0125_peaklist.h5')
    this_looper.get_tracks()

def check_particles(ds):
    bad_index=[]
    for grid in ds.index.grids:
        pos = grid['particle_position']
        for i in [0,1,2]:
            check=pos[:,i] > grid.RightEdge[i]
            check=np.logical_or(check,pos[:,i] < grid.LeftEdge[i])
            if check.any():
                locations=np.where(check)
                bad_index+= list( grid['particle_index'][check].v)
    return bad_index
#g=check_particles(this_looper.ds_list[120])
bad_list = []
for frame in list(range(122))+[125]:
    ds = yt.load("%s/DD%04d/data%04d"%(directory,frame,frame))
    this_bad_list = check_particles(ds)
    print("Frame %d len bad %d"%(frame,len(this_bad_list)))
    bad_list += this_bad_list
    bad_list = list(np.unique(bad_list))

