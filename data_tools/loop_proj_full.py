

from go import *

all_cores  =looper.get_all_nonzero()
all_cores.sort()
ni = int(sys.argv[1])
core_list= all_cores[ni:ni+10]
print("CORE_LIST", core_list)
print(len(core_list))

rm = rainbow_map(len(all_cores))
def cbump(value):
    output = (value+0).v
    output[output>1] = output[output>1]-1
    output[output<0] = output[output<0]+1
    return output

@looper.frame_loop
def full_proj(self, axis_list=[0,1,2]):
    for axis in axis_list:
        proj = self.ds.proj('density',axis,center='c')
        pw = proj.to_pw(center = 'c',width=(1.0,'code_length'), origin='domain')
        pw.set_cmap('density','gray')
        radius_from_core = []
        core_label = ""
        for nc,core_number in enumerate(self.target_indices):
            this_snapshot = self.make_snapshot(self.current_frame,core_number)
            if this_snapshot.R_centroid is None:
                this_snapshot.get_all_properties()
            ind =  np.where(core_number == all_cores)[0][0]
            color = rm(ind)
            print("COLOR", color,ind) 
            core_label += "c%04d_"%core_number
            pw.annotate_text(cbump(this_snapshot.R_centroid),
                             "%d"%core_number,text_args={'color':color}, 
                             inset_box_args={'visible':False},
                             coord_system='data')
            pw.annotate_select_particles(1.0, col=color, indices=self.target_indices[core_number])
        outname = '%s_full_particles_%sn%04d'%(self.out_prefix,core_label,self.current_frame)
        print( pw.save(outname))
looper.core_looper.full_proj = full_proj

if 1:
    directory = '/scratch2/dcollins/Paper19_48/B02/u05-r4-l4-128'
    this_looper = looper.core_looper(directory= directory,
                                     sim_name = 'u05',
                                     out_prefix = 'test',
                                     target_frame = 125,
                                     frame_list = [0,1,2]+list(range(10,130,5))+[125],
                                     core_list = core_list, #[THIS_CORE],
                                     fields_from_grid=['x'] #,'y','z']
                                  )
    this_looper.get_target_indices(h5_name='u05_0125_peaklist.h5',
                                     bad_particle_list='bad_particles.h5')
    this_looper.full_proj(axis_list=[0])
