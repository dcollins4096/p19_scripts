
from starter2 import *

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
#ni = int(sys.argv[1])
core_list=[79]# all_cores[ni:ni+10]
print("CORE_LIST", core_list)
print(len(core_list))
print(np.where(core_list==154))
rm = rainbow_map(len(all_cores))
#make an instance
#for THIS_CORE in core_list
if 1:
    directory = '/scratch2/dcollins/Paper19_48/B02/u05-r4-l4-128'
    this_looper = looper.core_looper(directory= directory,
                                     sim_name = 'u05',
                                     out_prefix = 'test',
                                     target_frame = 125,
                                     frame_list = list(range(11,120))+[125],
                                     core_list = core_list, #[THIS_CORE],
                                     fields_from_grid=['x'] #,'y','z']
                                  )
    this_looper.get_target_indices(h5_name='u05_0125_peaklist.h5',
                                     bad_particle_list='bad_particles.h5')
    #this_looper.get_tracks()

#core_31_density_baddies=nar([1499257, 1597829, 1618375, 1622656, 1626883, 1631042, 1634500])
#core_31_density_baddies=nar([1495097])
core_31_density_baddies=nar([1499257])
@looper.particle_loop
def select_particles(looper,snapshot,axis_list=[0,1,2]):
    for axis in axis_list:
        proj = snapshot.ds.proj('density',axis,center='c')
        pw = proj.to_pw(center = 'c',width=(1.0,'code_length'), origin='domain')
        pw.set_cmap('density','gray')
        radius_from_core = []
        core_label = ""
        these_particles = core_31_density_baddies
        pw.annotate_select_particles(1.0, col='r', indices=these_particles, p_size=5)
        pw.annotate_sphere(snapshot.R_centroid,snapshot.R_mag.max() ) #R_mag.max())
        outname = '%s_one_bad_%sn%04d'%(looper.out_prefix,snapshot.core_id,snapshot.frame)
        print(outname)
        print( pw.save(outname))
looper.core_looper.select_particles=select_particles
#this_looper.select_particles()

@looper.particle_loop   #the decorator takes care of the loop
def print_centroid(looper,snapshot): #here's a function, needs to take at least these arguments)
    print("Core %d frame %d centroid (%s)"%(
          snapshot.core_id, snapshot.frame, str( snapshot.R_centroid)))
looper.core_looper.print_centroid = print_centroid #add your function to looper.core_looper

#this_looper.print_centroid()

#this decrator makes the following function get called
#for all frames in core_looper.frame_list.


@looper.particle_loop
def core_proj_follow(looper,snapshot, field='density', axis_list=[0,1,2], color='r'):
    for ax in axis_list:
        scale_min = snapshot.ds.arr(0.05,'code_length')
        scale = max([snapshot.R_mag.max(), scale_min])
        sph = snapshot.ds.sphere(snapshot.R_centroid,scale)
        proj = snapshot.ds.proj(field,ax,center=snapshot.R_centroid, data_source = sph)
        pw = proj.to_pw(center = snapshot.R_centroid,width=(1.0,'code_length'), origin='domain')
        pw.zoom(1./scale.v)
        #pw = proj.to_pw(center = 'c',width=(1.0,'code_length'), origin='domain')
        pw.set_cmap(field,'gray')
        #pw.annotate_sphere(snapshot.R_centroid,snapshot.R_mag.max(), circle_args={'color':color} ) #R_mag.max())
        #pw.annotate_text(snapshot.R_centroid,
        #                 "%d"%snapshot.core_id,text_args={'color':color}, 
        #                 inset_box_args={'visible':False},
        #                 coord_system='data')
        pw.annotate_select_particles(1.0, col=color, indices=snapshot.target_indices)
        outname = "%s_c%04d_n%04d_centered"%(looper.out_prefix,snapshot.core_id,snapshot.frame)
        print(pw.save(outname))
    return pw
looper.core_looper.core_proj_follow = core_proj_follow
this_looper.core_proj_follow(axis_list=[0])
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

#this_looper.core_proj_follow()
#this_looper.full_proj()
#this_looper.core_test()
#this_looper.core_circle()
