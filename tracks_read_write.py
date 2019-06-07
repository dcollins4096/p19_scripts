
from starter1 import *
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
import time_kludge
reload(time_kludge)
#many1: frame_list = list(range(10,130,10)) + [125]
#many1: core_list =  [ 120, 10, 308, 41, 44, 110, 1, 65], 
if 0:
	if 'this_looper' not in dir():
		directory = '/scratch2/dcollins/Paper19_48/B02/u05-r4-l4-128'
		this_looper = looper.core_looper(directory= directory,
										 sim_name = 'u05',
										 out_prefix = 'test',
										 target_frame = 125, frame_list =  [1], #[0,1,2]+list(range(10,130,10))+[125],
										 core_list =  [10],
										 fields_from_grid=['x','y','z','velocity_magnitude','magnetic_field_strength',
														  'velocity_divergence'],
									  )
		this_looper.get_target_indices(h5_name='u05_0125_peaklist.h5')
		this_looper.get_tracks()
#this_looper.save_loop("test.h5")
	self = this_looper

def save_loop(self,fname):
    fptr = h5py.File(fname,'w')
    #just a bunch of values
    try:
        for value in [ "current_frame",
                       "data_template",
                       "sim_name",     
                       "frame_list",   
                       "core_list",    
                       "directory",    
                       "target_frame", 
                       "out_prefix",   
                      ]:
            result = self.__dict__[value]
            if result is None:
                continue
            fptr.create_dataset(value,data=result)
        fptr.create_dataset('fields_from_grid',data=np.array( [np.string_(s) for s in self.fields_from_grid]))
        ti_grp=fptr.create_group('target_indices')
        for core in self.target_indices:
            ti_grp.create_dataset(str(core),data=self.target_indices[core])

        fptr.create_group('snaps')
        for frame in self.snaps:
            frame_group_name = 'frame %d'%frame
            grp=fptr['snaps'].create_group(frame_group_name)
            for core in self.snaps[frame]:
                snp=grp.create_group("core %d"%core)
                this_snap = self.snaps[frame][core]
                for val in ["core_id",
                            "frame",
                            "mask",
                            "pos",
                            "R_centroid",
                            "R_vec",
                            "R_mag",
                            "N_vec",
                            "V_bulk",
                            "V_rel",
                            "V_radial"]:
                    result = this_snap.__dict__[val]
                    if result is not None:
                        snp.create_dataset(val,data=result)
                fv=snp.create_group('field_values')
                for val in this_snap.field_values:
                    fv.create_dataset(val, data=this_snap.field_values[val])

        tr_gr=fptr.create_group('track_manager')
        for val in ['particle_ids','core_ids','frames','times']:
            tr_gr.create_dataset(val,data= self.tr.__dict__[val])
        tr_di = tr_gr.create_group('track_dict')
        for k in self.tr.track_dict:
            tr_di[k] = self.tr.track_dict[k]
        #tr_gr[] = self.core_ids
        #tr_gr[] = self.frames
        #tr_gr[] = self.times
    except:
        raise
    finally:
        fptr.close()
def load_loop(self,fname):
    fptr = h5py.File(fname,'r')
    #just a bunch of values
    try:
        for value in [ "current_frame",
                       "data_template",
                       "sim_name",     
                       "frame_list",   
                       "core_list",    
                       "target_frame", 
                       "out_prefix",   
                       "fields_from_grid"
                      ]:
            if value in fptr:
                self.__dict__[value] = fptr[value].value
        if self.directory is None:
            self.directory = fptr['directory'].value
        #ti_grp=fptr.create_group('target_indices')
        #for core in self.target_indices:
        #    ti_grp.create_dataset(str(core),data=self.target_indices[core])
        for core in fptr['target_indices']:
            self.target_indices[int(core)]=fptr['target_indices'][core].value

        for frame_name in fptr['snaps']:
            frame_number = int(frame_name.split()[1])
            self.load(frame_number,dummy=True)

            for core_name in fptr['snaps'][frame_name]:
                core_number = int(core_name.split()[1])
                core_grp = fptr['snaps'][frame_name][core_name]
                self.snaps[frame_number][core_number]=self.make_snapshot(frame_number,core_number,
                                                                        dummy_ds=True)
                this_snap=self.snaps[frame_number][core_number]
                if 'time' in core_grp:
                    this_snap.time = core_grp['time'].value
                else:
                    this_snap.time = time_kludge.d[frame_number]
                #this_snap.field_values={}
                for val in ["core_id",
                            "frame",
                            "mask",
                            "N_vec"]:
                    if val in core_grp:
                        this_snap.__dict__[val]=core_grp[val].value
                for val in ["R_centroid",
                            "R_vec",
                            "pos",
                            "R_mag", ]:
                    if val in core_grp:
                        #thisval = this_snap.ds.arr(core_grp[val].value,'code_length')
                        thisval = yt.units.yt_array.YTArray(core_grp[val].value,'cm')
                        this_snap.__dict__[val]=thisval
                for val in ["V_bulk",
                            "V_rel",
                            "V_radial"]:
                    if val in core_grp:
                        #thisval = this_snap.ds.arr(core_grp[val].value,'code_velocity')
                        thisval = yt.units.yt_array.YTArray(core_grp[val].value,'cm/s')
                        this_snap.__dict__[val]=thisval
                for val in core_grp['field_values']:
                    this_snap.field_values[val]=core_grp['field_values'][val].value
            if self.tr is None:
                self.tr = trackage.track_manager(self)
            self.tr.ingest(this_snap)
        #for val in ['particle_ids','core_ids','frames','times']:
        #    self.tr.__dict__[val] =fptr['track_manager'][val].value
        #self.tr.track_dict={}
        #for k in fptr['track_manager']['track_dict']:
        #    self.tr.track_dict[k]=fptr['track_manager']['track_dict'][k].value

    except:
        raise
    finally:
        fptr.close()
#def check(one,two):
#	for fld in ['core_list',
#				'current_frame',
#				'data_template',
#				'directory',
#				'ds',
#			'ds_list',
#			'field_values',
#			'fields_from_grid',
#			'filename',
#			'frame_list',
#			'get_current_frame',
#			'get_region',
#			'get_target_indices',
#			'get_tracks',
#			'individual_particle_tracks',
#			'load',
#			'make_snapshot',
#			'out_prefix',
#			'sim_name',
#			'snaps',
#			'target_frame',
#			'target_indices',
#			'tr']:

#    this_looper.get_tracks()
#    this_looper.tr.write('cores_many_take3.h5')
def check(lpr):
    f0=lpr.frame_list[1]
    c0=lpr.core_list[1]
    for frame in lpr.frame_list[1:]:
        for core in lpr.core_list:
            if core not in lpr.snaps[frame]:
                print("missing d%d c%d"%(frame,core))
                continue
            for field in lpr.snaps[f0][c0].field_values:
                if field is 'V_radial':
                    continue
                fv0 = lpr.snaps[f0][core].field_values[field]
                #this_looper.snaps[0][77].field_values['density'].size
                fv1 = lpr.snaps[frame][core].field_values[field]
                if fv0.size != fv1.size:
                    print("n %d c %d f %s z %d %d"%(frame,core,field, fv0.size, fv1.size))
#check(this_looper)
                
        

#fname = "test2.h5"
#save_loop(self,fname)
#new = looper.core_looper()
#load_loop(new,fname)
#new.fields_from_grid = list(new.fields_from_grid)+['Bx']
#new.get_tracks()
#new.frame_list=[7]
#new.get_tracks()
