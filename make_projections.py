from starter2 import *
import xtra_energy
directory = '/scratch1/dcollins/Paper19/u17_sphere_linear'
this_looper = looper.core_looper(directory=directory,savefile="u17_test.h5")
if 'this_looper' not in dir() and False:

    directory = '/scratch2/dcollins/Paper19_48/B02/u05-r4-l4-128/GravPotential'
    this_looper = looper.core_looper(directory= directory,
                                     derived=[xtra_energy.add_test_energies],
                                     sim_name = 'u05',
                                     out_prefix = 'test',
                                     target_frame = 125,
                                     frame_list = [40],
                                     core_list = [79],
                                     fields_from_grid=['x','y','z']
                                  )
    this_looper.get_target_indices(h5_name='u05_0125_peaklist.h5',
                                     bad_particle_list='bad_particles.h5')
    #this_looper.get_tracks()

@looper.frame_loop
def proj_onecore(self, axis_list=[0,1,2],core_list=[], field='density'):
    for axis in axis_list:
        #center=self.ds.arr(nar([0.4,0.7,0.8]),'code_length')
        center = self.ds.arr(nar([0.5]*3),'code_length')
        #self.ds.periodicity = (True,True,True)
        self.proj = self.ds.proj(field,axis,center=center)
        self.proj = yt.SlicePlot(self.ds, axis=axis, fields=[field], center=center)


        #pw = self.proj.to_pw(center = center,width=(1.0,'code_length'), origin='domain')
        #pw.set_cmap(field,'gray')
        radius_from_core = []
        core_label = ""
        for nc,core_number in enumerate(core_list):
            core_label += "c%04d_"%core_number
            self.proj.annotate_select_particles(1.0, col='r', indices=self.target_indices[core_number])
            outname = '%s_full_particles_%sn%04d'%(self.out_prefix,core_label,self.current_frame)
        print( self.proj.save(outname))
looper.core_looper.proj_onecore=proj_onecore
#this_looper.frame_list=list(range(41,47))
proj_onecore(this_looper,axis_list=[0],core_list=[0],field='density')
