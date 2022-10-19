# ----------------------------------------------------------------------------------------
# 2D LASER-PLASMA
import math
from math import pi,sqrt
import numpy as np 
# CONSTANTS
c                   = 299792458.
lambda_SI           = 0.8e-6
omega_SI            = 2.*pi*c/lambda_SI
micron              = 1.e-6  * omega_SI/c
fs                  = 1.e-15 * omega_SI
l0                  = 2*pi  # laser wavelength
t0                  = l0    # optical cycle
				

dx = dy = micron/25.
Lx = Ly = 50.*micron 
npatchx = npatchy = 2
T_sim = 300.*fs

cfl = 0.98   
dt = cfl / np.sqrt((1./dx)**2+(1./dy)**2) 

every_optical_cycle = int(t0/dt)
every_fs = int(fs/dt)
N_steps = int(T_sim / dt) 

every_field = 10*every_fs 
every_track = 10*every_fs 
every_scalar = 5*every_fs 

fwhm = 30.*fs 
center = 50.*fs
a0 = 5.
waist = 3.*micron

n0e = 1e-2

A=1.
Z=1.


T = 0.1

start = 12*micron

nppc = 1

def ele_dens(x,y):
    if x>=start: return n0e
    else: return 0.
    
def ion_dens(x,y):
    if x>=start: return n0e/Z
    else: return 0.

#def my_filter(particles):
#    return (particles.y<0.5*Ly+waist)*(particles.y>0.5*Ly-waist) 
    
Main(
    geometry = "2Dcartesian",
    cell_length = [dx,dy],
    grid_length  = [Lx,Ly],
    timestep = dt,
    simulation_time = T_sim,
    interpolation_order = 2,
    EM_boundary_conditions = [ ['silver-muller'],['periodic']],    
    solve_poisson = False ,
    maxwell_solver = 'Yee',
    number_of_patches = [npatchx, npatchy],
    patch_arrangement = 'hilbertian',
    random_seed = smilei_mpi_rank,
    print_every = int(N_steps/100.0),
    reference_angular_frequency_SI = omega_SI
)

LaserGaussian2D(
    box_side         = "xmin",
    a0               = a0,
    omega            = 1.,
    focus            = [start, Ly*0.5],
    waist            = waist,
    incidence_angle  = 0.,
    polarization_phi = 0.,
    ellipticity      = 0.,
    time_envelope    = tgaussian(fwhm=fwhm, center=center)
)

Species(
    name      = "ele",
    position_initialization = "random",
    momentum_initialization = "maxwell-juettner",
    charge= -1.,
    particles_per_cell = nppc,
    mass = 1.,
    number_density = ele_dens,
    mean_velocity = [0.],
    temperature = [T],
    boundary_conditions = [['reflective', 'reflective'], ['periodic', 'periodic']],
    is_test = False,
    pusher = "boris",
)

Species(
    name      = "ion",
    position_initialization = "ele",
    momentum_initialization = "maxwell-juettner",
    charge= Z,
    particles_per_cell = nppc,
    mass = 1836.*A,
    number_density = ion_dens,
    mean_velocity = [0.],
    temperature = [T],
    boundary_conditions = [['reflective', 'reflective'], ['periodic', 'periodic']],
    is_test = False,
    pusher = "boris",
)



DiagScalar(every=every_scalar)

DiagFields(
    every = every_field,
    fields = ['Bz', 'Ex', 'Rho_ele', 'Rho_ion']
)

DiagTrackParticles(
    species = "ele",
    every = every_track,
    flush_every = 100,
    #filter = my_filter,
    attributes = ["x", "y", "px", "py", "pz", "w"] #, "Ex", "Ey", "Bz"]
)


#DiagTrackParticles(
#    species = "ion",
#    every = every_track,
#    flush_every = 100,
#    filter = my_filter,
#    attributes = ["x", "y", "px", "py", "pz", "w"] #, "Ex", "Ey", "Bz"]
#)

