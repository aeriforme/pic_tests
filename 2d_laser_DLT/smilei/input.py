import math
from math import pi,sqrt,sin,tan,log10,exp,ceil,log,atan
import numpy as np
from numpy import random, vectorize, trapz

#################
# PARAMETERS ####
#################

# GENERAL SI PARAMETERS
c = 299792458.
lambda_SI = 0.8e-6
omega_SI = 2.*pi*c/lambda_SI
micron = 1.e-6 * omega_SI/c
fs = 1.e-15 * omega_SI
l0 = 2.*pi # wavelength 
t0 = 2.*pi # optical cycle 
mc2 = 0.510998950e6 # MeV

# BOX 
dx = micron/20.
Lx = 70*micron 
nx = Lx/dx 
npatch_x = 8

dy = micron/20.
Ly = 30*micron 
ny = Ly/dy
npatch_y = 8

# TIME 
cfl = 0.98
dt = (cfl / sqrt( 1./dx**2 + 1./dy**2 )) 
T_sim = 1.5*Lx
#T_sim = 0.1*Lx

# DIAGNOSTICS
every_optical_cycle = int(t0/dt)
every_fs = int(fs/dt)
every_track = every_fs #only for test species 
N_steps = int(T_sim / dt) 

# LASER 
a0 = 20.
waist = 3.*micron
delay_peak = 60.*fs 
intensity_fwhm = 30.*fs
field_fwhm = intensity_fwhm*sqrt(2)
laser_length = 2*field_fwhm
x_spot = laser_length
y_spot = 0.5*Ly
 
# PLASMA
ne_f = 0.1
ne_s = 30.
ne_c= 50
foam_thick=20*micron;
subs_thick=1*micron;
cont_thick=0.1*micron;
x_foam_start = laser_length
x_foam_end = x_foam_start + foam_thick
x_subs_start = x_foam_end
x_subs_end = x_subs_start + subs_thick
x_cont_start = x_subs_end
x_cont_end = x_cont_start + cont_thick

#deuterium
Z_f=1.
A_f=2.
#boron-11
Z_s=5.
A_s=11.
#protons
Z_c=1.
A_c=1.

nppc = 8. #number of particles per cell

temp = 10./mc2


def ele_dens_foam(x,y) :
    if x < x_foam_start:
    	return 0.
    elif x >= x_foam_start and x < x_foam_end :
        return ne_f
    else : return 0.

def ele_dens_subs(x,y) :
    if x < x_subs_start:
    	return 0.
    elif x >= x_subs_start and x < x_subs_end :
        return ne_s
    else : return 0.
    
def ele_dens_cont(x,y) :
    if x < x_cont_start:
    	return 0.
    elif x >= x_cont_start and x < x_cont_end :
        return ne_c
    else : return 0.

def ion_dens_foam(x,y) :
    if x < x_foam_start:
    	return 0.
    elif x >= x_foam_start and x < x_foam_end :
        return ne_f/Z_f
    else : return 0.
    
def ion_dens_subs(x,y) :
    if x < x_subs_start:
    	return 0.
    elif x >= x_subs_start and x < x_subs_end :
        return ne_s/Z_s
    else : return 0.
        
def ion_dens_cont(x,y) :
    if x < x_cont_start:
    	return 0.
    elif x >= x_cont_start and x < x_cont_end :
        return ne_c/Z_c
    else : return 0.
            
       
    
###############
# SETUP #######
###############

Main(
    geometry = '2Dcartesian',
    interpolation_order = 2,
    number_of_patches = [npatch_x,npatch_y],
    simulation_time = T_sim,
    timestep = dt,
    reference_angular_frequency_SI = omega_SI ,
    cell_length = [dx,dy],
    grid_length = [Lx,Ly] ,
    EM_boundary_conditions = [['silver-muller'],['periodic']],
    print_every = int(T_sim/dt/100.),
    random_seed = smilei_mpi_rank,
    patch_arrangement = 'hilbertian',
    maxwell_solver = 'Yee',
    solve_poisson = False ,
)

# LOAD BALANCING
LoadBalancing(
    initial_balance = True,
    every = 0,
    cell_load = 1.,
    frozen_particle_load = 0.1
)

# VECTORIZATION
Vectorization(
    mode = 'off',
    reconfigure_every 	= 20,
)

LaserGaussian2D(
    box_side            = 'xmin',
    a0                  = a0 ,
    omega		= 1., 
    waist               = waist,
    focus               = [x_spot, y_spot],
    polarization_phi 	= 0.,
    ellipticity         = 0.,   # Linearly polarized
    incidence_angle     = 0. ,
    time_envelope       = tgaussian(fwhm=field_fwhm,center=delay_peak) 
)

Species(
    name = 'ele_f' ,
    position_initialization = 'random' ,
    momentum_initialization = 'maxwell-juettner' ,
    temperature = [temp],
    mean_velocity = [0.],
    particles_per_cell = nppc,
    mass = 1.,
    charge = -1.,
    number_density = ele_dens_foam,
    boundary_conditions =[['remove'],['periodic']],
    pusher = "boris"
)

Species(
    name = 'ele_f_test' ,
    position_initialization = 'random' ,
    momentum_initialization = 'maxwell-juettner' ,
    temperature = [temp],
    mean_velocity = [0.],
    particles_per_cell = 1,
    mass = 1.,
    charge = -1.,
    number_density = ele_dens_foam,
    boundary_conditions =[['remove'],['periodic']],
    is_test = True,
    pusher = "boris"

)

Species(
    name = 'ele_s' ,
    position_initialization = 'random' ,
    momentum_initialization = 'maxwell-juettner' ,
    temperature = [temp],
    mean_velocity = [0.],
    particles_per_cell = nppc,
    mass = 1.,
    charge = -1.,
    number_density = ele_dens_subs,
    boundary_conditions =[['remove'],['periodic']],
    pusher = "boris"
)


Species(
    name = 'ele_s_test' ,
    position_initialization = 'random' ,
    momentum_initialization = 'maxwell-juettner' ,
    temperature = [temp],
    mean_velocity = [0.],
    particles_per_cell = 1,
    mass = 1.,
    charge = -1.,
    number_density = ele_dens_subs,
    boundary_conditions =[['remove'],['periodic']],
    pusher = "boris",
    is_test = True
)


Species(
    name = 'ele_c' ,
    position_initialization = 'random' ,
    momentum_initialization = 'maxwell-juettner' ,
    temperature = [temp],
    mean_velocity = [0.],
    particles_per_cell = nppc,
    mass = 1.,
    charge = -1.,
    number_density = ele_dens_cont,
    boundary_conditions =[['remove'],['periodic']],
    pusher = "boris"
)


Species(
    name = 'ion_f',
    position_initialization = 'random' ,
    momentum_initialization = 'cold',
    mean_velocity = [0.],
    particles_per_cell = nppc ,
    mass = A_f * 1836. ,
    charge = Z_f ,
    number_density = ion_dens_foam,
    boundary_conditions =[['reflective'],['periodic']],
    pusher = "boris"
)

Species(
    name = 'ion_s',
    position_initialization = 'random' ,
    momentum_initialization = 'cold',
    mean_velocity = [0.],
    particles_per_cell = nppc ,
    mass = A_f * 1836. ,
    charge = Z_f ,
    number_density = ion_dens_subs,
    boundary_conditions =[['reflective'],['periodic']],
    pusher = "boris"
)


Species(
    name = 'ion_c',
    position_initialization = 'random' ,
    momentum_initialization = 'cold',
    mean_velocity = [0.],
    particles_per_cell = nppc ,
    mass = A_c * 1836. ,
    charge = Z_c ,
    number_density = ion_dens_cont,
    boundary_conditions =[['reflective'],['periodic']],
    pusher = "boris"
)

Species(
    name = 'ion_c_test',
    position_initialization = 'random' ,
    momentum_initialization = 'cold',
    mean_velocity = [0.],
    particles_per_cell = 1 ,
    mass = A_c * 1836. ,
    charge = Z_c ,
    number_density = ion_dens_cont,
    boundary_conditions =[['reflective'],['periodic']],
    pusher = "boris",
    is_test = True
)





###############
# OUTPUTS #####
###############

DiagScalar(every=every_fs)

DiagFields(
    every =10*every_fs,
    fields = ['Bz', 'Rho_ele_f','Rho_ele_s','Rho_ele_c','Rho_ion_f','Rho_ion_s','Rho_ion_c']
)


#DiagTrackParticles(
#    species = "ion_c_test",
#    every = every_track,
#    flush_every = 100,
#    #filter = my_filter,
#    attributes = ["x", "y", "px", "py", "pz", "w"] #, "Ex", "Ey", "Bz"]
#)


DiagParticleBinning(
    deposited_quantity = "weight",
    every = every_fs,
    time_average = 1,
    species = ["ele_f","ele_s","ele_c"],
    axes = [ ["px",   -1,    1,    100] ]
)

DiagParticleBinning(
    deposited_quantity = "weight",
    every = every_fs,
    time_average = 1,
    species = ["ion_f","ion_s","ion_c"],
    axes = [ ["px",   -10,    10,    100] ]
)


DiagParticleBinning(
    deposited_quantity = "weight",
    every = 5,
    time_average = 1,
    species = ["ele_f","ele_s","ele_c","ion_f","ion_s","ion_c"],
    axes = [ ["ekin",    0.1*mc2,	40*mc2,   100, "logscale"] ]
)

