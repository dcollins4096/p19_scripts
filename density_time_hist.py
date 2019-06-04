import matplotlib.colors as colors
if 0:
    nbins=100
    thhist = np.zeros([nbins,len(asort)])
    mind = density.min()
    maxd = density.max()
    bins = np.logspace(np.log10(mind), np.log10(maxd),nbins+1)
    fig,ax = plt.subplots(1,1)
    for itime,ntime in enumerate(asort):
        thhist[:,itime]=np.histogram(density[:,ntime],bins)[0]
    norm = colors.LogNorm(vmin=1,vmax=thhist.max())
    ploot=ax.imshow(thhist,origin='lower',interpolation='nearest',norm=norm)
    ax.set_aspect(.3)
    cbar = plt.colorbar(ploot,norm=norm)
    cbar.set_clim(vmin=1,vmax=thhist.max())
    cbar.cmap.set_under('w')

    fig.savefig('density_time_hist.png')
if 1:
    fig2,ax2=plt.subplots(1,1)
    nparticles = density.shape[0]
    times_ordered=thtr.times[asort]
    for nt,t in enumerate(times_ordered):
        ax2.scatter([t]*nparticles,density[:,asort[nt]],c='k',s=.01)
    ax2.set_yscale('log')
    ax2.set_ylim(0.01,1e8)
    fig2.savefig('density_time_scatter.png')
        
plt.close('all')
