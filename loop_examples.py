
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

directory = '/scratch1/dcollins/Paper19/SphereTest/s09_periodic_sphere'

#this decrator makes the following function get called
#for all frames in core_looper.frame_list.
@looper.frame_loop
def full_proj(self, axis_list=[0,1,2]):
    for axis in axis_list:
        proj = self.ds.proj('density',axis,center='c')
        pw = proj.to_pw(center = 'c',width=(1.0,'code_length'), origin='domain')
        pw.set_cmap('density','gray')
        radius_from_core = []
        core_label = ""
        for nc,core_number in enumerate(self.target_indices):
            core_label += "c%04d_"%core_number
            pw.annotate_select_particles(1.0, col='r', indices=self.target_indices[core_number])
        outname = '%s_full_particles_%sn%04d'%(self.out_prefix,core_label,self.current_frame)
        print( pw.save(outname))
looper.core_looper.full_proj = full_proj


@looper.particle_loop
def core_proj_follow(looper,snapshot, axis_list=[0,1,2], color='r'):
    for ax in axis_list:
        proj = snapshot.ds.proj('density',ax,center=snapshot.R_centroid)
        pw = proj.to_pw(center = snapshot.R_centroid,width=(1.0,'code_length'), origin='domain')
        pw.set_cmap('density','gray')
        pw.annotate_sphere(snapshot.R_centroid,snapshot.R_mag.max(), circle_args={'color':color} ) #R_mag.max())
        xax = looper.ds.coordinates.x_axis[ax]
        yax = looper.ds.coordinates.x_axis[ax]
        pw.annotate_text(snapshot.R_centroid,
                         "%d"%snapshot.core_id,text_args={'color':color}, 
                         inset_box_args={'visible':False},
                         coord_system='data')
        pw.annotate_select_particles(1.0, col=color, indices=snapshot.target_indices)
        outname = "%s_c%04d_n%04d_centered"%(looper.out_prefix,snapshot.core_id,snapshot.frame)
        print(pw.save(outname))
    return pw
looper.core_looper.core_proj_follow = core_proj_follow

@looper.particle_loop
def core_circle(self, axis_list=[0,1,2]):
    for axis in axis_list:
        proj = self.ds.proj('density',axis,center='c')
        pw = proj.to_pw(center = 'c',width=(1.0,'code_length'), origin='domain')
        pw.set_cmap('density','gray')
        radius_from_core = []
        core_label = ""
        pw.annotate_select_particles(1.0, col='r', indices=self.target_indices[self.core_id])
        pw.annotate_sphere(self.R_centroid,self.R_mas.max() ) #R_mag.max())
        outname = '%s_fullcircle_cores_%sn%04d'%(self.out_prefix,self.core_id,self.current_frame)
        print( pw.save(outname))
looper.core_circle=core_circle


@looper.particle_loop
def core_test(looper,snapshot):
    snapshot.do_the_stuff()
    print("FRAME %d CORE_ID %d N_PARICLES %s Centroid %s"%(
        snapshot.frame, snapshot.frame, str(snapshot.pos.shape), str(snapshot.R_centroid)))
looper.core_looper.core_test = core_test
#a

fname = 'run2.py'
exec(open(fname).read(), globals(), locals())

#this_looper.full_proj()
#this_looper.core_test()
#this_looper.core_circle()
