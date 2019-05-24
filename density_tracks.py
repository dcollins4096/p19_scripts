
import matplotlib
matplotlib.use('Agg')
import yt
import math
import matplotlib.pyplot as plt
import numpy as np
import pdb
nar = np.array
from collections import defaultdict

from importlib import reload

import looper
reload(looper)
import trackage
reload(trackage)
import tracks_read_write as trw
reload(trw)
from davetools import *
plt.close('all')
file_list=glob.glob('/home/dcollins/scratch/Paper19/particle_error/particle_error_test_c0031_threeframes.h5')
file_list=glob.glob('/home/dcollins/scratch/Paper19/particle_error/track_indfix_sixteenframe_core_0031.h5')
file_list=glob.glob('/scratch1/dcollins/Paper19/Datasets/track_indfix_sixteenframe_core_*.h5')
directory = '/home/dcollins/scratch/u05-r4-l4-128'
class parameter_package():
    def __init__(self,rho0=1,rho1=1,tc=1,tff_local=1,rho_c=1, tff_global=1):
        self.rho0   =rho0
        self.rho1   =rho1
        self.tc     =tc
        self.tff_local=tff_local
        self.tff_global=tff_global
        self.rho_c  =rho_c
    def set(self,rho0=1,rho1=1,tc=1,tff_local=1,rho_c=1, tff_global=1):
        self.rho0   =rho0
        self.rho1   =rho1
        self.tc     =tc
        self.tff_global=tff_global
        self.tff_local=tff_local
        self.rho_c  =rho_c
    def __str__(self):
        output = ""
        output += "rho0     =%0.2e\n"%rho0
        output += "rho1     =%0.2e\n"%rho1
        output += "tc       =%0.2e\n"%tc
        output += "tff_local=%0.2e\n"%tff_local
        output += "tc/tffg  =%0.2e\n"%(tc/tff_global)
        output += "rho_c    =%0.2e\n"%rho_c
        return output
freefall = defaultdict(parameter_package)

for nfile,fname in enumerate(file_list[:3]) :#[:3])
    #0164.h5
    t1 = fname.split("/")[-1]
    #l = len("track_three_to_test_core_")
    #l = len("track_sixteen_frames_core_")
    l = len("track_indfix_sixteenframe_core_")

    this_cor = int(t1[l:l+4]) #[fname.index('_'):]
    #this_cor=31
    #if this_cor not in  [12, 31]:
    #    continue
    print(this_cor)
    this_looper=looper.core_looper(directory=directory)
    trw.load_loop(this_looper,fname)
    thtr = this_looper.tr
    all_cores = np.unique(thtr.core_ids)
    core_list=all_cores
    rm = rainbow_map(len(all_cores))

    if 1:
        #time plots
        asort =  np.argsort(thtr.times)
        tsorted = thtr.times[asort]
        for nc,core_id in enumerate(core_list):
            ms = trackage.mini_scrubber(thtr,core_id)
            density = thtr.c([core_id],'density')
            tmap=rainbow_map(ms.ntimes)
            plt.clf()
            for npart in range(ms.nparticles):
                plt.plot( tsorted, density[npart,asort],c='k',linestyle=':',marker='*')


            t0 = thtr.times[asort][0]
            t1 = thtr.times[asort][-1]
            rho0 = np.mean(density[:,asort[0]])
            rho1 = np.mean(density[:,asort[-1]])
            alpha = 1.8
            tc = t1*(1-(rho1/rho0)**(-1./alpha))**-0.5
            G=1620./(4*np.pi)
            tff_global = np.sqrt(3*np.pi/(32*G*1))
            tff_local = np.sqrt(3*np.pi/(32*G*rho0))
            rhot = rho0*(1-(tsorted/tc)**2)**-alpha
            rho_c = 3*np.pi/(32*G*tc**2)
            freefall[core_id].set(rho0=rho0,rho1=rho1,tc=tc, tff_global=tff_global,
                                  tff_local=tff_local,rho_c=rho_c)
            plt.plot( tsorted, rhot, c='r')
            plt.text( tsorted[0], rho1, r'$tc = %0.2e \rho_c = %0.2e$'%(tc,rho_c))
            plt.text( tsorted[0], 0.5*rho1, r'$tc/tff = %0.2e$'%(tc/tff_global))
            outname = 'image_tracks/rho_t_fit_c%04d.png'%(core_id)
            plt.yscale('log')
            plt.savefig(outname)
            print('saved '+outname)

