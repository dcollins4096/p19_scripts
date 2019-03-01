
from go import *
# for new ds.all_data() profile plot - 
import core_dump
def rs(arr):
    return arr.reshape(arr.shape + (1,)) 
def rs2(arr):
    return arr.reshape(arr.shape + (1,) + (1,)) 

import trackage
reload(trackage)
if 'track' not in dir():
    track = trackage.track_manager(None)

    directory = '/scratch2/dcollins/Paper19_48/B02/u05-r4-l4-128'
    ds0 = yt.load(directory+"/DD0000/data0000")
    ad0 = ds0.all_data()
    extrema = {'density':[5e-3,100]}
    n_bins=128

    #track.read('all_large_0000_0125.h5')
    #track.read('tracks_first_last_all_nonzero.h5')
    track.read('tracks_frame012_all_nonzero.h5')
all_cores = np.unique(track.core_ids)
if 1:
    plt.clf()
    for i,nc in enumerate(all_cores):
        this_density = track.c(int(nc),'velocity_magnitude')[:,0]
        #plt.hist(np.log10(this_density),histtype='step',color='k')
        plt.hist(this_density,histtype='step',color='k')
    plt.yscale('log')
    #plt.xlabel(r'$\log_{10} ||v||$')
    plt.xlabel(r'$||v||$')
    plt.ylabel(r'$N(v)$')
    plt.savefig('/home/dcollins4096/PigPen/p56_allpdf_vel_notlog.png')
if 0:
    plt.clf()
    for i,nc in enumerate(all_cores):
        this_density = track.c(int(nc),'density')[:,0]
        plt.hist(np.log10(this_density),histtype='step',color='k')
    plt.yscale('log')
    plt.savefig('/home/dcollins4096/PigPen/p47_allhist.png')


if 0:
    means = np.zeros_like(all_cores,dtype='float')
    stds = np.zeros_like(all_cores,dtype='float')
    npart = np.zeros_like(all_cores,dtype='float')
    for i,nc in enumerate(all_cores):
        this_density = np.log(track.c(int(nc),'density')[:,0])
        means[i]=this_density.mean()
        npart[i] = this_density.size
        stds[i]=np.sqrt(np.sum((this_density-means[i])**2)/npart[i])
    plt.clf()
    plt.scatter(means,stds)
    plt.ylabel(r'$\sqrt{\langle \ln \rho^2 \rangle}$')
    plt.xlabel(r'$\langle \ln \rho \rangle$')
    plt.savefig('/home/dcollins4096/PigPen/p47_log_mean_std.png')

