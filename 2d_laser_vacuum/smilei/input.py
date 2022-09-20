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
mc2 = 0.510998950e6 # eV

# BOX 
dx = micron/20.
Lx = 102.4*micron 
nx = Lx/dx 
npatch_x = 64

dy = micron/20.
Ly = 51.2*micron 
ny = Ly/dy
npatch_y = 64

# TIME 
cfl = 0.98
dt = (cfl / sqrt( 1./dx**2 + 1./dy**2 )) 
T_sim = 1.5*Lx

# DIAGNOSTICS
every_optical_cycle = int(t0/dt)
every_fs = int(fs/dt)
N_steps = int(T_sim / dt) 
 

# LASER 
a0 = 20.
waist = 3.*micron
delay_peak = 60.*fs 
intensity_fwhm = 30.*fs
field_fwhm = intensity_fwhm*sqrt(2)
x_spot = 0.5*Lx 
y_spot = 0.5*Ly

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
    EM_boundary_conditions = [['PML'],['PML']],
    number_of_pml_cells = [[10,10],[10,10]],
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


###############
# OUTPUTS #####
###############

DiagScalar(every=2*every_fs)

DiagFields(
    every =10*every_fs,
    fields = ['Ey']
)
