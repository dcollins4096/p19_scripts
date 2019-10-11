from starter2 import *
import matplotlib.image as mpimg

import data_locations as dl
reload(dl)

#file_list=glob.glob('%s/track_indfix_sixteenframe_core_*.h5'%dl.snapshot_location)
#file_list=glob.glob('/home/dcollins/scratch/Paper19/all_frames/track_all_frames_*.h5')
gloob='/home/dcollins/scratch/Paper19/all_frames/track_all_frames_c0079.h5'
gloob='/Volumes/deep7/RESEARCH2/Paper19_47_overlap/2019-07-24-sphere/*h5'

file_list=glob.glob(gloob)
plt.close('all')
output_dir='sphere'


for nfile,fname in enumerate(file_list) :#[:3])

    #this is crude
    #t1 = fname.split("/")[-1]
    #l = len("track_indfix_sixteenframe_core_")
    #this_cor = int(t1[l:l+4]) #[fname.index('_'):]
    #if this_cor not in  [12]:#, 31]:
    #    continue
    #print(this_cor)

    this_looper=looper.core_looper(directory=dl.enzo_directory)
    trw.load_loop(this_looper,fname)
    thtr = this_looper.tr
    all_cores = np.unique(thtr.core_ids)
    core_list=all_cores
    rm = rainbow_map(len(all_cores))

    if 1:
        #time plots
        asort =  np.argsort(thtr.times)
        n0=asort[0]
        tsorted = thtr.times[asort]
        fig,ax=plt.subplots(1,1)
        for nc,core_id in enumerate(core_list):
            ms = trackage.mini_scrubber(thtr,core_id)
            density = thtr.c([core_id],'density')
            tmap=rainbow_map(ms.ntimes)
            norm = mpl.colors.Normalize()
            norm.autoscale( np.log10(density[:,n0]))
            cmap = mpl.cm.jet
            color_map = mpl.cm.ScalarMappable(norm=norm,cmap=cmap)

            for npart in range(ms.nparticles):# [::100]:
                c = color_map.to_rgba(np.log10(density[npart,n0]))
                #ax.plot( tsorted, density[npart,asort],c='k',linestyle=':',marker='*')
                ax.plot( tsorted, density[npart,asort],c=c,linewidth=.1)#linestyle=':')
            #err= np.exp(np.log(density).std(axis=0)[asort])
            ax.plot(tsorted, density.mean(axis=0)[asort],c='k')
            #ax.errorbar(tsorted, density.mean(axis=0)[asort],c='k',yerr=err)
            #ax.plot(tsorted, density.mean(axis=0),c='k')

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
            ax.plot( tsorted, rhot, c='r',label=r'$tc/tff = %0.2e$'%(tc/tff_global))

            ax.plot( tsorted, rhot, c='r',label=r'$tc/tff = %0.2e$'%(tc/tff_global))
            #ax.text( tsorted[0], rho1, r'$tc = %0.2e \rho_c = %0.2e$'%(tc,rho_c))
            #ax.text( tsorted[0], 0.5*rho1, 
            ax.legend(loc=0)
            axbonk(ax,xscale='linear',yscale='log',xlabel='t',ylabel=r'$\rho$')
            name_base = "%s/rho_time/"%output_dir
            oname = name_base + "/density_time_c%04d"%(core_id)
            fig.savefig(oname)
            print('saved '+oname)
#           for i,n in enumerate(asort):
#               timeline=plt.plot( [tsorted[i]]*2,[1,1e8],c=[0.5]*3,linewidth=0.1)
#               timetext=plt.text( tsorted[i], 1e8, 'n=%d'%thtr.frames[n])
#               outname = name_base+'/rho_t_fit2_c%04d_n%04d.png'%(core_id,i)
#               plt.yscale('log')
#               plt.savefig(outname)
#               timeline[0].remove()
#               timetext.remove()

#               print('saved '+outname)

#               frame = thtr.frames[n]
#               frame = i #you suck.
#               basedir = "/home/dcollins/RESEARCH2/Paper19_47_overlap/0000_density_tracks/x/ALL"
#               core_dir = "%s/c%04d"%(basedir,core_id)
#               image = "%s/test_full_particles_c%04d_n%04d_Projection_x_density.png"%(core_dir,core_id,frame)
#               img1 = mpimg.imread(image) 
#               img2 = mpimg.imread(outname) 
#               both = np.zeros( [ max([img1.shape[0],img2.shape[0]]), img1.shape[1]+img2.shape[1],4])
#               both[ 0:img1.shape[0] , 0:img1.shape[1]]=img1
#               both[ 0:img2.shape[0] , img1.shape[1]:]=img2
#               oname2 = "image_tracks/density_2_c%04d_n%04d"%(core_id,frame)
#               print(oname2)
#               plt.imsave(oname2,both)

