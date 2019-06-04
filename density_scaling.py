from go import *
import davetools
reload(davetools)
extents = davetools.extents
file_list=glob.glob('/home/dcollins/scratch/Paper19/track_index_fix/track_indfix_sixteenframe_core_*.h5')
plt.close('all')
    
if 'rho_extents' not in dir():
    rho_extents=extents()
    x_extents=extents()
    r_extents=extents()
    alpha_extents=extents()
    for nfile,fname in enumerate(file_list):
        this_looper=looper.core_looper(directory=directory)
        trw.load_loop(this_looper,fname)
        all_cores = np.unique(thtr.core_ids)
        for nc,core_id in enumerate(core_list):
            ms = trackage.mini_scrubber(thtr,core_id)
            if ms.nparticles == 1:
                continue
            density = thtr.c([core_id],'density')
            t2 = np.tile(thtr.times,(ms.nparticles,1))
            alpha = density*t2**2
            r = ms.r
            x = r[t2>0]/t2[t2>0]
            rho_extents(density)
            alpha_extents(alpha[t2>0])
            r_extents(r)
            x_extents(x)


for nfile,fname in enumerate(file_list[:2]):
    #0164.h5
    t1 = fname.split("/")[-1]
    #l = len("track_three_to_test_core_")
    #l = len("track_sixteen_frames_core_")
    l = len("track_indfix_sixteenframe_core_")

    this_cor = int(t1[l:l+4]) #[fname.index('_'):]
    #this_cor=31
    #if this_cor not in  [12]:#, 31]:
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
        n0=asort[0]
        tsorted = thtr.times[asort]
        fig=plt.figure(figsize=(8,4))
        axa=fig.subplots(1,2)
        #fig,axa=plt.subplots(1,2, figsize=(8,8))
        ax=axa[1]
        ax2=axa[0]
        #fig2,ax2=plt.subplots(1,1)

        for nc,core_id in enumerate(core_list):
            ms = trackage.mini_scrubber(thtr,core_id)
            if ms.nparticles == 1:
                continue
            density = thtr.c([core_id],'density')
            tmap=rainbow_map(ms.ntimes)
            norm = mpl.colors.Normalize()
            norm.autoscale( np.log10(density[:,n0]))
            cmap = mpl.cm.jet
            for n_count,n_time in enumerate(asort):
                time=thtr.times[n_time]
                c=tmap(n_count,ms.nparticles)
                this_r=ms.r[:,n_time]
                r_un = nar(sorted(np.unique(this_r)))
                ax2.scatter(this_r,density[:,n_time],c=c,label=thtr.times[n_time],s=0.1)
                ax2.plot(r_un, 100*(r_un/1e-2)**-2,c='k',linewidth=0.1)

                this_x = this_r#/thtr.times[n_time]
                alpha = density[:,n_time]*this_r**2#/time**2
                ax.scatter(this_x,alpha,c=c,label=thtr.times[n_time],s=0.1)
                x_un=nar(sorted(np.unique(this_x)))
                ax.plot(x_un, 1e-2*x_un**0,c='k',linewidth=0.1)
            davetools.axbonk(ax,xscale='log',yscale='log',xlabel='x',ylabel=r'$\alpha$')
                             #xlim=x_extents.minmax, ylim=alpha_extents.minmax)
            davetools.axbonk(ax2,xscale='log',yscale='log',xlabel='r',ylabel=r'$\rho$',
                             xlim=r_extents.minmax, ylim=rho_extents.minmax)
            outname = 'image_tracks/density_rescale_c%04d'%core_id
            fig.savefig(outname)
            #outname = 'image_tracks/density_scale_c%04d'%core_id
            #fig2.savefig(outname)

            print("saved",outname)
