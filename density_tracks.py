from starter2 import *

import data_locations as dl
reload(dl)

file_list=glob.glob('%s/track_indfix_sixteenframe_core_*.h5'%dl.snapshot_location)
plt.close('all')


for nfile,fname in enumerate(file_list) :#[:3])

    #this is crude
    t1 = fname.split("/")[-1]
    l = len("track_indfix_sixteenframe_core_")
    this_cor = int(t1[l:l+4]) #[fname.index('_'):]
    if this_cor not in  [12]:#, 31]:
        continue
    print(this_cor)

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
        for nc,core_id in enumerate(core_list):
            ms = trackage.mini_scrubber(thtr,core_id)
            density = thtr.c([core_id],'density')
            tmap=rainbow_map(ms.ntimes)
            plt.clf()
            norm = mpl.colors.Normalize()
            norm.autoscale( np.log10(density[:,n0]))
            cmap = mpl.cm.jet
            color_map = mpl.cm.ScalarMappable(norm=norm,cmap=cmap)

            for npart in range(ms.nparticles):
                c = color_map.to_rgba(np.log10(density[npart,n0]))
                #plt.plot( tsorted, density[npart,asort],c='k',linestyle=':',marker='*')
                plt.plot( tsorted, density[npart,asort],c=c,linewidth=.1)#linestyle=':')
            err= np.exp(np.log(density).std(axis=0)[asort])
            #plt.plot(tsorted, density.mean(axis=0)[asort],c='k')
            plt.errorbar(tsorted, density.mean(axis=0)[asort],c='k',yerr=err)
            #plt.plot(tsorted, density.mean(axis=0),c='k')




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
            plt.plot( tsorted, rhot, c='r')
            plt.text( tsorted[0], rho1, r'$tc = %0.2e \rho_c = %0.2e$'%(tc,rho_c))
            plt.text( tsorted[0], 0.5*rho1, r'$tc/tff = %0.2e$'%(tc/tff_global))
            for i,n in enumerate(asort[0:1]):
                timeline=plt.plot( [tsorted[i]]*2,[1,1e8],c=[0.5]*3,linewidth=0.1)
                timetext=plt.text( tsorted[i], 1e8, 'n=%d'%thtr.frames[n])
                outname = 'image_tracks/rho_t_fit2_c%04d_s%04d.png'%(core_id,i)
                plt.yscale('log')
                plt.savefig(outname)
                timeline[0].remove()
                timetext.remove()

                print('saved '+outname)

