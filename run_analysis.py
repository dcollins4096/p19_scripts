
import matplotlib
matplotlib.use('Agg')
import yt
import math
import matplotlib.pyplot as plt
import numpy as np
nar = np.array
import time
import pyximport; pyximport.install()
import particle_grid_mask
#from p19b_select_particle_callback import *
from p19_tools import *
from yt.analysis_modules.level_sets.api import *
import copy 
import pdb


"""
This is a re-factored script for analyzing tracer particles in clumps.
Structure:
    1.) This routine is the parent routine that 
    """

#Get indices
if 'leaf_indices' not in dir():
    if 1:
        #this top section needs to define
        #sim, scratchdir, data_template, leaf_indices.
        sim = 'u05'
        late_frame = 125
        scratchdir ='/scratch2/dcollins/Paper19_48/B02/u05-r4-l4-128'
        fname = '%s/DD%04d/data%04d'%(scratchdir,late_frame,late_frame)
        late_ds = yt.load(fname)
        keepers = [0,1,8,10,11,12,67,68,64,61, 201, 125, 306]
        keepers = [34]
        keepers=[10]
        print("READ")
        leaf_indices = get_leaf_indices(late_ds,pickle_name = 'u05_0125_peaklist.pickle', 
                                     subset = keepers)
        data_template = '%s/DD%04d/data%04d'
        #loc_2 = nar([ 0.45874023,  0.48266602,  0.20336914] )
        #peaks = fPickle.load('u05_0125_peaklist.pickle')
    if 0:
        #this top section needs to define
        #sim, scratchdir, data_template, leaf_indices.
        sim = 's09'
        scratchdir = '/scratch1/dcollins/Paper19/SphereTest/s09_periodic_sphere'
        data_template = '%s/DD%04d/data%04d'
        first_ds = yt.load(data_template%(scratchdir,0,0))
        first_sphere = first_ds.sphere([0.5]*3, 0.25)
        leaf_indices = {0:first_sphere['particle_index']}


if 'framelist' not in dir():
    framelist = [2]
    individual_particle_tracks = False
    stop = False
    test_output = False
    do_proj = False
    analysis_loop=True
    do_field_stuff=True
    fields_from_grid = ['density', 'cell_volume','x','y','z', 'velocity_x', 'velocity_y','velocity_z','dx']
    other_args={}
    field_values_by_frame={}

import stuff_to_do
reload(stuff_to_do)
do_stuff = stuff_to_do.radial_velocity

if test_output:
    output=do_stuff(ds,sim,frame,ic, pos,vel,field_values, stop, other_args)
if analysis_loop:
    for frame in framelist:
        print("New Frame",frame)
        fname = data_template%(scratchdir,frame,frame)
        ds = yt.load(fname)
        G_code = ds['GravitationalConstant']
        ad = ds.all_data()
        #times.append(ds['InitialTime'])
        #cycles.append(ds['InitialCycleNumber'])

        #projection
        if do_proj:
            for axis in [0,1,2]:
                proj = ds.proj('density',axis,center='c')
                pw = proj.to_pw(center = 'c',width=(1.0,'code_length'), origin='domain')
                pw.set_cmap('density','gray')
                radius_from_clump = []
                clump_label = ""
                for nc,clump_number in enumerate(leaf_indices):
                    clump_label += "c%04d_"%clump_number
                    pw.annotate_select_particles(1.0, col='r', indices=leaf_indices[clump_number])
                outname = '%s_with_clumps_%sn%04d'%(sim,clump_label,frame)
                print pw.save(outname)

    
        if do_field_stuff:
            for ic in leaf_indices: #you can loop over dictionaries.
                found_any, mask = get_indices(leaf_indices[ic], ad)
                pos = ad['particle_position'][mask == 1]
                vel = ad['particle_velocity'][mask == 1]
                ind = ad['particle_index'][mask == 1]
                field_values={}
                for field in fields_from_grid:
                    field_values[field] = np.zeros(leaf_indices[ic].size)-1
                get_particle_values_from_grid(ds,field_values, leaf_indices[ic], pos)

                if  (field_values['cell_volume'] < 0).any():
                    NM= (field_values['cell_volume'] < 0).sum()
                    print("ERROR: some particles (%d of them) not found.  This is problematic."%NM)

                plt.clf()
                output=do_stuff(ds,sim,frame,ic, leaf_indices, pos,vel,field_values, stop, other_args)




