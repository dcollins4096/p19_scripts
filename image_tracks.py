
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
import trackage
reload(trackage)
import tracks_read_write as trw
reload(trw)
from davetools import *
if 0:
    #newer read.
    if 'this_looper' not in dir():
        directory = '/home/dcollins/scratch/u05-r4-l4-128'
        this_looper=looper.core_looper(directory=directory)
        #file_list=glob.glob('/home/dcollins/scratch/Paper19/track_sixteen/*h5')
        file_list=glob.glob('/home/dcollins/scratch/Paper19/track_three/*h5')
        #file_list=glob.glob('/home/dcollins/scratch/Paper19/track_sixteen_good/*h5')

        for fname in file_list:

            trw.load_loop(this_looper,fname)
        thtr = this_looper.tr
        
import copy
def shift_up(pos,point):
    out = np.sort(copy.copy(pos),axis=1)
    out[ out < point ] = out[ out < point]+1
    return out
bork=0
delta=0
bo=0
smo=0
def shift_down(pos):
    #shift based on the time history: if a particle jumps more than half the box,
    #move it.
    global bork
    global delta
    global bo
    global smo
    sign = -1
    out = copy.copy(pos)
    delta = sign*(out[:,1:]-out[:,:-1])
    shape = delta.shape
    bork=delta.max(axis=1)
    bork = np.tile(bork,(shape[1],1)).transpose()
    bo = np.logical_and(delta <= bork, delta > 0.5)

    out[:,:-1][bo] -= 1
   

    distance_from_final = np.abs(out- np.tile(out[:,-1], (out.shape[1],1)).transpose())
    ft = np.abs(distance_from_final[:,:-1]) > 0.5
    out[:,:-1][ft] += 1*np.sign(distance_from_final[:,:-1][ft])



    #out[ out > point ] = out[ out > point]-1
    return out#,delta

if 0:
    #first read 
    if 'this_looper' not in dir():
        #directory = '/scratch2/dcollins/Paper19_48/B02/u05-r4-l4-128'
        directory = '/home/dcollins/scratch/u05-r4-l4-128'
        this_looper = looper.core_looper(directory= directory,
                                         sim_name = 'u05',
                                         out_prefix = 'test',
                                         target_frame = 125,
                                         frame_list = list(range(10,130,10)) + [125],
                                         core_list =  [ 120, 10, 308, 41, 44, 110, 1, 65], 
                                         fields_from_grid=['x','y','z']
                                      )
    if this_looper.tr is None:
        this_looper.tr = trackage.track_manager(this_looper)
        this_looper.tr.read('cores_many1.h5')
        thtr=this_looper.tr

#this_raw_x=raw_x[549,:]
#this_raw_y=raw_y[549,:]
#b=shift_down(this_raw_x)
##b=shift_down(this_y[0:1,:])
    ls=all_cores
#ls = [202]
#ls=[41]
#ls = [10]
file_list=glob.glob('/home/dcollins/scratch/Paper19/track_three/*h5')
directory = '/home/dcollins/scratch/u05-r4-l4-128'
for fname in file_list[:3]:
    #0164.h5
    t1 = fname.split("/")[-1]
    l = len("track_three_to_test_core_")
    this_cor = int(t1[l:l+4]) #[fname.index('_'):]
    print(this_cor)
    if this_cor not in [25]:
        continue
    this_looper=looper.core_looper(directory=directory)
    trw.load_loop(this_looper,fname)
    thtr = this_looper.tr
    all_cores = np.unique(thtr.core_ids)
    ls=all_cores
    rm = rainbow_map(len(all_cores))
    plt.clf()
    plt.plot([0,1,1,0,0],[0,0,1,1,0])
    if 1:
        #just make images.
        fig_rhist,ax_rhist=plt.subplots(1,1)
        fig,ax=plt.subplots(1,1,figsize=(8,8))
        fig_rmst,ax_rmst=plt.subplots(1,1,figsize=(8,8))
        ax.plot([0,1,1,0,0],[0,0,1,1,0])
        for nc,core_id in enumerate(ls):
            ax.clear()
            print('plot core %d'%core_id)
            n_for_color = int(np.where( all_cores == core_id)[0])
            raw_x = thtr.c([core_id],'y')
            raw_y = thtr.c([core_id],'z')
            if 1:
                #do the shift
                this_x = shift_down(raw_x)
                this_y = shift_down(raw_y)
            else:
                #don't actuall shift
                this_x = raw_x+0
                this_y = raw_y+0
            mean_x = np.mean(this_x,axis=0)
            mean_y = np.mean(this_y,axis=0)
            nparticles,ntimes=this_x.shape
            meanx2 = np.tile(mean_x,(raw_x.shape[0],1))
            meany2 = np.tile(mean_y,(raw_x.shape[0],1))
            r2 = (this_x-meanx2)**2 + (this_y-meany2)**2
            r=np.sqrt(r2)
            rmax = np.max(r,axis=0)
            max_track = np.where( r[:,0] == rmax[0])
            rmax_fat=np.tile(rmax,(raw_x.shape[0],1))
            rms = np.sqrt( np.mean(r2,axis=0))
            ax_rmst.plot(thtr.times, rms,c=rm(n_for_color))
            ok = np.where( r==rmax_fat)
            circle_max = plt.Circle( (mean_x[0],mean_y[0]), rmax[0], color=rm(n_for_color),fill=False)
            ax.add_artist(circle_max)
            circle_rms = plt.Circle( (mean_x[0],mean_y[0]), rms[0], color=rm(n_for_color),fill=False)
            ax.add_artist(circle_rms)
            #theseparts = np.arange(0,nparticles,100,dtype='int')
            #theseparts = [413]
            theseparts = np.arange(0,nparticles,dtype='int')
            rr=rainbow_map(len(theseparts))
            ax.plot([0,1,1,0,0],[0,0,1,1,0])
            for npp,npart in enumerate(theseparts):
                lab=None
                if npart==0:
                    lab = 'c %d'%core_id
                    #print('wtf',lab)
                    this_color = rm(n_for_color)
                if 0:
                    this_color = rr(npp)
                ax.plot( this_x[npart,:], this_y[npart,:],c=this_color,label=lab)
            delta = 0.5
            ax.set_xlim(-delta,1+delta); ax.set_ylim(-delta,1+delta)
            fig.savefig('image_tracks/image_c%04d.png'%core_id)

            ntimes = len(raw_x[0,:]) 
            tmap=rainbow_map(ntimes)
            ax_rhist.clear()
            for nt in range(ntimes):
                ax_rhist.hist(r[:,nt],histtype='step',color=tmap(nt))
            ax_rhist.set_title('core %d'%core_id)
            fig_rhist.savefig('image_tracks/rhist_c%04d.png'%core_id)

            #lab = 'c %d'%core_id
            ##ax.plot(mean_x , mean_y,c=rm(nc),label=lab)
            #ax.plot(this_x-meanx2 ,this_y-meany2,c=rm(nc),label=lab)

        ax.set_xlim(-1,2); ax.set_ylim(-1,2)
        #ax.set_xlim(-delta,1+delta); ax.set_ylim(-delta,1+delta)
        ax.legend(loc=0)
        fig.savefig('image_tracks/core_rel.png')
        fig_rmst.savefig('image_tracks/rrms.png')
        plt.close('all')

    if 1:
        #radial plots.
        #fig_rhist,ax_rhist=plt.subplots(1,1)
        fig,ax=plt.subplots(1,1,figsize=(8,8))
        #fig_rmst,ax_rmst=plt.subplots(1,1,figsize=(8,8))
        ax.plot([0,1,1,0,0],[0,0,1,1,0])
        for nc,core_id in enumerate(ls):
            n_for_color = int(np.where( all_cores == core_id)[0])
            raw_x = thtr.c([core_id],'y')
            raw_y = thtr.c([core_id],'z')
            density = thtr.c([core_id],'density')
            this_x = shift_down(raw_x)
            this_y = shift_down(raw_y)
            mean_x = np.mean(this_x,axis=0)
            mean_y = np.mean(this_y,axis=0)
            nparticles,ntimes=this_x.shape
            meanx2 = np.tile(mean_x,(raw_x.shape[0],1))
            meany2 = np.tile(mean_y,(raw_x.shape[0],1))
            r2 = (this_x-meanx2)**2 + (this_y-meany2)**2
            r=np.sqrt(r2)
            ntimes = len(raw_x[0,:]) 
            tmap=rainbow_map(ntimes)
            if 0:
                #density plots
                ntimes = len(raw_x[0,:]) 
                tmap=rainbow_map(ntimes)
                plt.clf()
                for nt in range(ntimes):
                    plt.scatter(r[:,nt], density[:,nt],c=[tmap(nt)]*r.shape[0])
                plt.xscale('log');plt.yscale('log')
                plt.xlim(1e-5,0.5)
                plt.ylim(0.01,4e6)
                plt.savefig('image_tracks/rho_t_c%04d.png'%core_id)
            if 0:
                #density plots
                ntimes = len(raw_x[0,:]) 
                tmap=rainbow_map(ntimes)
                plt.clf()
                for nt in range(ntimes):
                    plt.scatter(r[:,nt], density[:,nt],c=[tmap(nt)]*r.shape[0])
                plt.xscale('log');plt.yscale('log')
                plt.xlim(1e-5,0.5)
                plt.ylim(0.01,4e6)
                plt.savefig('image_tracks/rho_t_c%04d.png'%core_id)



