from starter2 import *
import xtra_operators as xo
reload(xo)
kinetic_validators=[yt.ValidateSpatial(1,['x-velocity','y-velocity','z-velocity','kinetic_energy'])]
std_validators = [yt.ValidateSpatial(1,['Bx','By','Bz', 'x-velocity','y-velocity','z-velocity'])]
std_validators_2 = [yt.ValidateSpatial(1,['magnetic_field_x','magnetic_field_y','magnetic_field_z', 'x-velocity','y-velocity','z-velocity'])]
kinetic_validators=[yt.ValidateSpatial(1,['x-velocity','y-velocity','z-velocity','kinetic_energy'])]
kinetic_validators=[yt.ValidateSpatial(1,['density'])]
pressure_validators=[yt.ValidateSpatial(1,['x-velocity','y-velocity','z-velocity','pressure'])]
def grad(data,fieldname,direction):
    iM1 = slice(None,-2)
    iP1 = slice(2,None)
    all = slice(1,-1)
    all_all=tuple([all]*3)
    dxi=1./(2*data.dds )
    out = np.zeros_like(data[fieldname]*dxi[0])
    Left = [all]*3
    Right = [all]*3
    Right[direction] = iP1
    Left[direction] = iM1
    Left=tuple(Left); Right=tuple(Right)
    out[all_all] = (data[fieldname][ Right ]- data[fieldname][ Left]) *dxi[direction]
    return out

def add_force_terms(obj):
    def momentum_flux(field,data):
        f1 =xo.gradf(0.5*data['x-velocity']*data['kinetic_energy'],0,data.dds)
        f2 =xo.gradf(0.5*data['y-velocity']*data['kinetic_energy'],1,data.dds)
        f3 =xo.gradf(0.5*data['z-velocity']*data['kinetic_energy'],2,data.dds)
        return (f1+f2+f3)
    obj.add_field('momentum_flux',momentum_flux, validators=kinetic_validators, sampling_type='cell',
                  units='dyne/(cm**2*s)')
def _grav_pot_grad(field,data):
    gx = grad(data,'PotentialField',0)
    gy = grad(data,'PotentialField',1)
    gz = grad(data,'PotentialField',2)
def add_energies(obj):
    def grav_energy(field,data):
        gx =data.ds.arr(grad(data,'PotentialField',0),'code_length/code_time**2')
        gy =data.ds.arr(grad(data,'PotentialField',1),'code_length/code_time**2')
        gz =data.ds.arr(grad(data,'PotentialField',2),'code_length/code_time**2')
        alpha = 1./data.ds['GravitationalConstant'] #=1/4 pi G
        alpha = data.ds.quan(alpha, '1/(code_length**3/code_time**2/code_mass)')
        return ( (gx**2+gy**2+gz**2)*alpha )
    obj.add_field('grav_energy',grav_energy,validators=[yt.ValidateGridType()],
                 units='code_mass*code_length**2/(code_time**2*code_length**3)')
    def therm_energy(field,data):
        sound_speed = data.ds.quan(1.,'code_velocity')
        rho_0 = data.ds.quan(1.,'code_mass/code_length**2')
        e = sound_speed**2*np.log( data['density']/rho_0)
        therme = data['density']*e
        return therme
    obj.add_field('therm_energy',therm_energy,
                 units='code_mass*code_length**2/(code_time**2*code_length**3)')
def add_test_energies(obj):
    def gx(field,data):
        gx =data.ds.arr(grad(data,'PotentialField',0),'code_length/code_time**2')
        return gx
    obj.add_field('gx',gx,validators=[yt.ValidateGridType()],
                 units='code_length/code_time**2')

def add_force_terms(obj):
    def grav_x(field,data):
        gi =-data.ds.arr(grad(data,'PotentialField',0),'code_length/code_time**2')
        return gi
    obj.add_field('grav_x',grav_x,validators=[yt.ValidateSpatial(1,['PotentialField'])],
                 units='code_length/(code_time**2)',sampling_type='cell')
    def grav_y(field,data):
        gi =-data.ds.arr(grad(data,'PotentialField',1),'code_length/code_time**2')
        return gi
    obj.add_field('grav_y',grav_y,validators=[yt.ValidateSpatial(1,['PotentialField'])],
                 units='code_length/(code_time**2)',sampling_type='cell')
    def grav_z(field,data):
        gi =-data.ds.arr(grad(data,'PotentialField',2),'code_length/code_time**2')
        return gi
    obj.add_field('grav_z',grav_z,validators=[yt.ValidateSpatial(1,['PotentialField'])],
                 units='code_length/(code_time**2)',sampling_type='cell')
    def pressure_work(field,data):
        try:
            gi =xo.AdotDel(data, ['velocity_x','velocity_y','velocity_z'], 'gas_pressure') #-data.ds.arr(grad(data,'PotentialField',2),'code_length/code_time**2')
        except:
            return np.zeros_like(data['density'])
        
        return gi
    work_units='dyne/(cm**2*s)'
    obj.add_field('pressure_work',pressure_work,validators=pressure_validators,
                 units=work_units,sampling_type='cell', take_log=True)
    def gravity_work(field,data):
        try:
            gi =data['x-velocity']*data['grav_x']+
                data['x-velocity']*data['grav_x']+
                data['x-velocity']*data['grav_x']);
        except:
            return np.zeros_like(data['density'])
        
        return gi
    work_units='dyne/(cm**2*s)'
    obj.add_field('gravity_work',gravity_work,validators=pressure_validators,
                 units=work_units,sampling_type='cell', take_log=True)
    def gas_pressure(field,data):
        return data['density']*data.ds.quan(1,'code_velocity')**2
    obj.add_field('gas_pressure',gas_pressure,units='dyne/cm**2')
    def momentum_flux(field,data):
        f1 = xo.gradf(0.5*data['x-velocity']*data['kinetic_energy'],0,data.dds)
        f2 = xo.gradf(0.5*data['y-velocity']*data['kinetic_energy'],1,data.dds)
        f3 = xo.gradf(0.5*data['z-velocity']*data['kinetic_energy'],2,data.dds)
        return (f1+f2+f3)
    obj.add_field('momentum_flux',momentum_flux, validators=kinetic_validators, sampling_type='cell',
                  units='dyne/(cm**2*s)')
    #def pressure_force(field,data):
    #    output  =-data['x-velocity']*xo.gradf(data['gas_pressure'],0,data.dds)
    #    output +=-data['y-velocity']*xo.gradf(data['gas_pressure'],1,data.dds)
    #    output +=-data['z-velocity']*xo.gradf(data['gas_pressure'],2,data.dds)
    #    output = data.ds.arr(np.zeros(data['density'].shape), 'dyne/(cm**2*s)')
    #    return output
    #obj.add_field('pressure_force',pressure_force,validators=pressure_validators,
    #              units='dyne/(cm**2*s)')

