from go import *
import os
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
#many1: frame_list = list(range(10,130,10)) + [125]
#many1: core_list =  [ 120, 10, 308, 41, 44, 110, 1, 65], 
#frame 120 core 96
core_list = all_nonzero.astype('int')[::-1]
#core_list=[31]
#if 'this_looper' not in dir():
for core in core_list:
    directory = '/scratch2/dcollins/Paper19_48/B02/u05-r4-l4-128'
    output_name = 'track_indfix_sixteenframe_core_%04d.h5'%core
    if os.path.exists(output_name):
        continue
    this_looper = looper.core_looper(directory= directory,
                                     sim_name = 'u05',
                                     out_prefix = 'test',
                                     target_frame = 125,
                                     frame_list = [0,1,2]+list(range(10,130,10))+[125],
                                     core_list =  [core],# core_list,
                                     fields_from_grid=['x','y','z','velocity_magnitude','magnetic_field_strength',
                                                      'velocity_divergence']
                                  )
    this_looper.get_target_indices(h5_name='u05_0125_peaklist.h5',
                                     bad_particle_list='bad_particles.h5')
    this_looper.get_tracks()
    trw.save_loop(this_looper,output_name)
    #this_looper.tr.write('big_partial_frame_0120_0125.h5')
#sn = this_looper.snaps[120][96]
"""Note on actualizing this data extraction:
    1.) It's huge and expensive, need to be able to write data to disk.
    2.) Need to make sure we're not making globally conserved quantities:
        multiply sampled cells!!! Need to paint an array and mask.
    3.) There are some problems with the shift operator and/or field values.


Individaul core issues:
    1.) frame 120 core 96:
        * Field values has 609 entries,
        * R_centroid has 607.
        * target_indices has 609.
        * Arg.  The particles in question are in fact off the edge of the grid that contains them.
    """

