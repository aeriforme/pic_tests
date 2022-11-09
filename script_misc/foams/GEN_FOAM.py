#!/bin/python 
import math
from math import pi,sqrt,sin,tan,log10,exp
import numpy as np
from numpy import random, vectorize, trapz
import matplotlib.pyplot as plt
from scipy.constants import * 
import sys

# GENERAL SI PARAMETERS AND l0 DEFINITION
lambda_SI = 0.8 * micron # wavelength
omega_SI = 2.0 * pi * c / lambda_SI
n_crit = 1.1 * 1e21 / centi**(3) * (1.e-6/lambda_SI)**(2) #m^-3

# BOX
Ly = 60*micron
Lx = 70.*micron 

# X
dx = micron / 62.5 # resolution
nx = int(Lx / dx) # Number of cells in X
# Y
dy = dx
ny = int(Ly / dy) # Number of cells in Y



# FOAM
foam_thick = Lfx = float(sys.argv[1])*micron
ave_ele_dens=float(sys.argv[2])*n_crit

Z_foam = 6 
ff = 0.025
np_ele_dens = ave_ele_dens/ff
foam_size_y = Lfy = Ly
sub_start_x = 0.7*Lx 
Fsx = sub_start_x - foam_thick
Fsy = -Ly*0.5

name = 'SPHERES_t'+str(int(foam_thick/micron))+'_d'+str(int(ave_ele_dens/n_crit))

d_x = 0.08*micron # diameter of the nanospheres
d_y = d_x
r = d_x/2. # radius of the nanospheres


Z = np.zeros([nx,ny]) # Matrix to store and show results

centres =  [ [ 0. , 0.] ] # List which will contain the coordinates of all the centres of the nanospheres
ff_flag = True
nanop_tot_volume = pi * r**2
foam_volume = Lfx * Lfy
foam_total_cells = foam_volume / (dx*dy)
target_ff = ff # FILLING FACTOR we want. will be used later to check we did a good job
cpnp = int(r**2 *pi / (dx*dy)) # Number of cells filled by each nanosphere (Cells Per NanoParticle)
filled_cells = 0.# Total number of cells filled with solid density material
failed_attempts = 0 # Just a variable to check the efficinecy of the script
while ff_flag : # If the filling factor is not reached it will keep adding nanospheres in the box
    next_x = random.random() * Lfx + Fsx # X cordinate of the next nanosphere
    next_y = random.random() * Lfy + Fsy # Y cordinate of the next nanosphere
    centres.append([next_x,next_y]) # Growing the list of the centres
    print len(centres) # Updating how many nanospheres we have in total
    success = True
    for c in range(0,len(centres)-1) : # Loop controlling that there is no overlapping among all already placed nanospheres and the newly added one
        if (centres[c][0] - centres[len(centres)-1][0]-dx)**2 + (centres[c][1] - centres[len(centres)-1][1]-dy)**2 <= r**2:
            del centres[-1] # If there is any overlapping delete the last nanosphere placed 
            success = False
            failed_attempts += 1
            break
    if success : # If there is no overlapping fill directly the matrix
        X = np.linspace(centres[-1][0]-r-2*dx,centres[-1][0]+r+2*dx,int(d_x/dx)+4) # We do not want to travel all the matrix everytime, but only
        Y = np.linspace(centres[-1][1]-r-2*dx,centres[-1][1]+r+2*dx,int(d_y/dx)+4) # the cells nearby the centre
        iz = int(X[0]/dx) # Index relative to the matrix "Z"
        jz = int(Y[0]/dx) # Index relative to the matrix "Z"
        for i in range(0, len(X)) :
            for j in range(0, len(Y)) : # Run through the matrix only around the last centre created
                if (Y[j]-centres[-1][1])**2 <= (r**2) * ( 1 - ((X[i]-centres[-1][0])**2)/r**2)  and iz+i < nx and jz+j <ny and Z[iz+i][jz+j]==0:
                    Z[iz+i][jz+j] = 1 # Fill the cells in the matrix
                    filled_cells += 1 # Update the count on the total cells filled
        ff = filled_cells / foam_total_cells # Update the filling factor
    if ff > target_ff : # If the filling factor is reached, then exit the cicle
        ff_flag = False

plt.imshow(Z,origin='lower') # Show the matrix
plt.colorbar()
plt.show(block=False)
np.save(name+'.npy',Z) # Save the matrix
plt.savefig(name+'.png',dpi=300)

f = open(name+'_ne.in', 'w+b')
byte_arr = np.transpose(Z)*np_ele_dens
binary_format = bytearray(byte_arr)
f.write(binary_format)
f.close()

f = open(name+'_ni.in', 'w+b')
byte_arr = np.transpose(Z)*np_ele_dens/Z_foam
binary_format = bytearray(byte_arr)
f.write(binary_format)
f.close()



################### CHECK ################################ This will count from zero all the filled cells to check consistency
filled_cells    = 0.
empty_cells     = 0.
total_cells     = 0.
ff_check        = 0.
for i in range(0,nx):
    for j in range(0,ny):
        print i+1,'/',nx
        if Z[i][j] == 1 :
            filled_cells += 1
        elif Z[i][j] == 0 :
            empty_cells += 1
ff_check = filled_cells/foam_total_cells
print 'ACTUAL FILLING FACTOR:\t',ff_check
print 'RELATIVE ERROR:\t',(ff_check-target_ff)/target_ff * 100,'%'
print 'FILLING FACTOR AFTER THE LOOP:\t', ff



print(nx, ny) 




