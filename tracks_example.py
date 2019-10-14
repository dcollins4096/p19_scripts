
import matplotlib
matplotlib.use('Agg')
import yt
import math
import matplotlib.pyplot as plt
import numpy as np
import pdb
nar = np.array

from importlib import reload

import looper
reload(looper)
if 'this_looper' not in dir():
    directory = '/scratch2/dcollins/Paper19_48/B02/u05-r4-l4-128'
    this_looper = looper.core_looper(directory= directory,
                                     sim_name = 'u05',
                                     out_prefix = 'test',
                                     target_frame = 125,
                                     frame_list = [0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,105,110,115,120,125],# list(range(10,130,10)) + [125],
                                     core_list =  [65], #, 167], 167 is broken...
                                     fields_from_grid=['x','y','z']
                                  )
    this_looper.get_target_indices(h5_name='u05_0125_peaklist.h5')
    this_looper.get_tracks()
    #this_looper.tr.write('test_file.h5')
#this_looper.tr.read('test_file.h5')
#
# the track manager is new, and stored in 
# looper.tr
# It takes all particles and makes time-series for them.
# looper.tr can be accessed like a dict:
# >>> density_tracks = looper.tr['density']
# and will return an array for the given field for each time.
# The array `density_tracks` has particles in the rows, and 
# time in the column.
# Other arrays of use are:
# >>> looper.tr.times
# >>> looper.tr.frames
# >>> looper.tr.core_ids
# >>> looper.tr.particle_ids

# >>> core_track = looper.tr.c(core_id,field)
# will return just the tracks for a field of  the particles
# in an individual core.

# >>> particle_track = looper.tr.p(particle_id,field)
# will produce tracks from a list of particles.

if 1:
    plt.clf()
    for track in this_looper.tr['density']:
        plt.plot(this_looper.tr.times, track, marker = '*')
    plt.xlabel('t')
    plt.ylabel(r'$\rho$')
    plt.yscale('log')
    plt.savefig('test2/many_tracks_many_cores.png')

for core_id in [25,65,100]: #this_looper.core_list:
    ds1 = yt.load(fname)
    ad = ds1.all_data()
    print(ds1['InitialTime'])
    G_code = 6.67E-8
    rho_0 = (ad['density'].v.sum())/len(ad['density'].v)
    tff = ((3*math.pi)/(32*G_code*rho_0))**0.5

    plt.clf()
    this_core_field = this_looper.tr.c(core_id,'density')
    for track in this_core_field:
        plt.plot(this_looper.tr.times, track, marker = '*')
    #Then you can make the mean and standard deviations
    mean = this_core_field.mean(axis=0)
    std = this_core_field.std(axis=0)
    plt.plot(this_looper.tr.times,mean, c='k',marker = '*')
#    plt.plot(this_looper.tr.times,mean+std, c='k',marker = '*')
#    plt.plot(this_looper.tr.times,mean-std, c='k',marker = '*')
    plt.plot(this_looper.tr.times,rho*(1-(this_looper.tr.times/tff)**2)**(-1.86),c = 'r')

    plt.xlabel('t')
    plt.ylabel(r'$\rho$')
    plt.title('Core %d'%core_id)
    plt.yscale('log')
    plt.savefig('test2/many_tracks_c%04d.png'%core_id)
