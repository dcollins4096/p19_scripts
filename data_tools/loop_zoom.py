
from go import *

# Do this with python 3.
# execfile is gone, do this instead:
# 
# >>> fname = 'thing.py'
# >>> exec(open(fname).read(), globals(), locals())

# The looper package has two main objects:
# looper.core_looper keeps track of the data, frames, target indices, frames, cores to examine.
# looper.snapshot hangs on to particle data for a given target core_id and frame.  
# There are two decorators, looper.frame_loop and looper.paticle_loop.
# The decorators take a function you write, and takes care of calling it for
# every core_id and frame in your instance of looper.core_looper.
# 
# So, if I want to print the centroid of each core, I'd do
all_cores  =looper.get_all_nonzero()
all_cores.sort()
core_list = [8]
print("CORE_LIST", core_list)
print(len(core_list))
print(np.where(core_list==154))
rm = rainbow_map(len(all_cores))
standard_sixteen = [0,1,2]+list(range(10,130,5))#+[125]
frame_list = standard_sixteen[-2:]
#make an instance
#for THIS_CORE in core_list
if 1:
    directory = '/scratch2/dcollins/Paper19_48/B02/u05-r4-l4-128'
    this_looper = looper.core_looper(directory= directory,
                                     sim_name = 'u05',
                                     out_prefix = 'zoom',
                                     target_frame = 125,
                                     frame_list = frame_list,
                                     core_list = core_list, #[THIS_CORE],
                                     fields_from_grid=['x'] #,'y','z']
                                  )
    this_looper.get_target_indices(h5_name='u05_0125_peaklist.h5',
                                     bad_particle_list='bad_particles.h5')

@looper.particle_loop
def core_proj_follow(looper,snapshot, axis_list=[0,1,2], color='r'):
    for ax in axis_list:
        proj = snapshot.ds.proj('density',ax,center=snapshot.R_centroid)
        pw = proj.to_pw(center = snapshot.R_centroid,width=(1.0,'code_length'), origin='domain')
        pw.zoom(2)
        pw.set_cmap('density','gray')
        pw.annotate_sphere(snapshot.R_centroid,snapshot.R_mag.max(), circle_args={'color':color} ) #R_mag.max())
        pw.annotate_text(snapshot.R_centroid,
                         "%d"%snapshot.core_id,text_args={'color':color}, 
                         inset_box_args={'visible':False},
                         coord_system='data')
        pw.annotate_select_particles(1.0, col=color, indices=snapshot.target_indices)
        print("wtf mate", snapshot.target_indices)
        outname = "%s_c%04d_n%04d_centered"%(looper.out_prefix,snapshot.core_id,snapshot.frame)
        print(pw.save(outname))
    return pw
looper.core_looper.core_proj_follow = core_proj_follow
this_looper.core_proj_follow(axis_list=[0])

