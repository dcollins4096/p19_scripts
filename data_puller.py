from starter2 import *
import data_locations as dl
all_nonzero = looper.get_all_nonzero()
#many1: frame_list = list(range(10,130,10)) + [125]
#many1: core_list =  [ 120, 10, 308, 41, 44, 110, 1, 65], 
#frame 120 core 96
core_list = all_nonzero.astype('int')[-3:]
frame_list = [1,125]
fields = ['x','y','z','velocity_magnitude','magnetic_field_strength', 'velocity_divergence']
output_base = "some_track_set"
#core_list=[31]
for core in core_list:
    output_name = '%s_c%04d.h5'%(output_base,core)
    if os.path.exists(output_name):
        continue
    this_looper = looper.core_looper(directory= dl.directory,
                                     sim_name = 'u05',
                                     out_prefix = 'test',
                                     target_frame = 125,
                                     frame_list = frame_list,
                                     core_list =  [core],# core_list,
                                     fields_from_grid=fields,
                                  )
    this_looper.get_target_indices(h5_name='u05_0125_peaklist.h5',
                                     bad_particle_list='bad_particles.h5')
    this_looper.get_tracks()
    trw.save_loop(this_looper,output_name)
