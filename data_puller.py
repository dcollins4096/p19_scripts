from starter2 import *
import data_locations as dl
all_nonzero = looper.get_all_nonzero()
#many1: frame_list = list(range(10,130,10)) + [125]
#many1: core_list =  [ 120, 10, 308, 41, 44, 110, 1, 65], 
#frame 120 core 96
if 0:
    core_list = all_nonzero.astype('int')
    frame_list = [0,1]+list(range(10,130,10))+[125]
    fields = ['x','y','z','velocity_magnitude','magnetic_field_strength', 'velocity_divergence']
    fields += ['velocity_x','velocity_y','velocity_z']
    fields += ['magnetic_field_%s'%s for s in 'xyz']
    output_base = "all_primitives"
    derived=[]
if 1:
    import xtra_energy
    core_list = all_nonzero.astype('int')
    frame_list = [0,1]+list(range(10,130,50))+[125]
    fields = ['mag_work']
    derived=[xtra_energy.add_force_terms]
    output_base = 'mag_work_tmp'
for core in core_list:  
    print("core %d"%core)
    output_name = '%s_c%04d.h5'%(output_base,core)
    if os.path.exists(output_name):
        print("File exists, skipping "+output_name)
        continue
    this_looper = looper.core_looper(directory= dl.enzo_directory,
                                     sim_name = 'u05',
                                     out_prefix = 'test',
                                     target_frame = 125,
                                     frame_list = frame_list,
                                     core_list =  [core],# core_list,
                                     fields_from_grid=fields,
                                     derived = derived
                                  )
    this_looper.get_target_indices(h5_name='u05_0125_peaklist.h5',
                                     bad_particle_list='bad_particles.h5')
    this_looper.get_tracks()
    trw.save_loop(this_looper,output_name)
