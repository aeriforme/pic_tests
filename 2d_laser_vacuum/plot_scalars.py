import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import micron,c,pi,centi,femto
import openpmd_api as io
import happi 
import sdf 
from matplotlib import use

use('AGG')

#__________________________________________________________
# create plot directory
plot_dir = './plots'
if os.path.exists(plot_dir) is False:
    os.mkdir(plot_dir)


#__________________________________________________________
lambda_SI = 0.8*micron # wavelength
omega_SI = 2.0*pi*c / lambda_SI
n_crit = 1.1*1e21 / centi**(3) * (1.e-6/lambda_SI)**(2) #m^-3
my_dpi = 300 

#__________________________________________________________
# simulation parameters
Lx = 102.4*micron 
resx = 20. # points per micron 
nx = Lx * resx / micron 
dx = Lx/nx

cfl = 0.98
dt = cfl * dx / (c * np.sqrt(2.))

Tsim = 1.5*Lx / c 
nsteps = int(Tsim/dt) 


every_fs = np.floor(femto/dt)
out_freq = int(10*every_fs)

steps = out_freq * np.arange(nsteps)
times = steps * dt/femto

steps = np.arange(0,nsteps,out_freq)

# warpx
warpx_dir='./warpx/diags/reducedfiles'
data_w = np.loadtxt(warpx_dir+'/FieldEnergy.txt')  
times_w = data_w[:,1]/femto
Uelm_w = data_w[:,2]
UE_w = data_w[:,3]
UB_w = data_w[:,4]

#smilei 
smilei_dir = './smilei'
s = happi.Open(smilei_dir) 

data_s = s.Scalar(scalar='Uelm', units=['J/m', 'fs'])
times_s = data_s.getTimes()
Uelm_s = data_s.getData()

UEx_s = s.Scalar(scalar='Uelm_Ex', units=['J/m', 'fs']).getData()
UEy_s = s.Scalar(scalar='Uelm_Ey', units=['J/m', 'fs']).getData()
UEz_s = s.Scalar(scalar='Uelm_Ez', units=['J/m', 'fs']).getData()

UE_s = np.asarray(UEx_s) + np.asarray(UEy_s) + np.asarray(UEz_s)

print(type(UE_s), type(UEx_s), np.shape(UE_s), np.shape(UEx_s))
UBx_s = s.Scalar(scalar='Uelm_Bx_m', units=['J/m', 'fs']).getData()
UBy_s = s.Scalar(scalar='Uelm_By_m', units=['J/m', 'fs']).getData()
UBz_s = s.Scalar(scalar='Uelm_Bz_m', units=['J/m', 'fs']).getData()
UB_s = np.asarray(UBx_s) + np.asarray(UBy_s) + np.asarray(UBz_s) 


# epoch
epoch_dir = './epoch/Data/'
Uelm_e=[]
times_e=[]

#sh.list_variables(data)
filenames=np.genfromtxt(epoch_dir+'scalars.visit',dtype='str')
ndumps=len(filenames) 

for t in range(ndumps):
    fname=epoch_dir+'scalars%04d.sdf' % t
    data=sdf.read(fname)
    times_e = np.append(times_e,data.Header['time']/femto)
    Uelm_e=np.append(Uelm_e,data.Total_Field_Energy_in_Simulation__J_.data)     


#__________________________________________________________
# plot images
my_dpi = 300 
fig, ax = plt.subplots(ncols=1, nrows=1, figsize=(1000./my_dpi, 1000./my_dpi), dpi=my_dpi)

ax.plot(times_w, Uelm_w, label='WarpX')
ax.plot(times_s, Uelm_s, label='Smilei')
ax.plot(times_e, Uelm_e, label='EPOCH')

ax.set_xlabel('time [fs]')
ax.set_ylabel('energy [J/m]')
ax.legend()
ax.set_title('field energy')

image_file_name =plot_dir+'/scalars.png' 
plt.tight_layout()
plt.savefig(image_file_name,dpi=my_dpi)
plt.close()


fig, ax = plt.subplots(ncols=1, nrows=1, figsize=(1000./my_dpi, 1000./my_dpi), dpi=my_dpi)

ax.plot(times_w, UE_w, label='E WarpX')
ax.plot(times_s, UE_s, label='E Smilei')
ax.plot(times_w, UB_w, label='B WarpX')
ax.plot(times_s, UB_s, label='B Smilei')

ax.set_xlabel('time [fs]')
ax.set_ylabel('energy [J/m]')
ax.legend()
ax.set_title('field energy')

image_file_name =plot_dir+'/scalars2.png' 
plt.tight_layout()
plt.savefig(image_file_name,dpi=my_dpi)
plt.close()



