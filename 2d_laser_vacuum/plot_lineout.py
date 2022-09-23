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
Lx = 70*micron 
resx = 25. # points per micron 
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
warpx_dir='./warpx/diags/Fields'
series = io.Series(warpx_dir+"/openpmd_%T.bp",io.Access.read_only)
#for n,ts in enumerate(series.iterations):

#smilei 
smilei_dir = './smilei'
s = happi.Open(smilei_dir) 
# Ey =  s.Field(0,'Ey', units=["um", "V/m", "fs"])
# for n, ts in Ey.getAvailableTimesteps():

# epoch
epoch_dir = './epoch/Data/'
filenames=np.genfromtxt(epoch_dir+'fields.visit',dtype='str') 
#for t in range(len(filenames)): 

for n,ts in enumerate(steps):
    # warpx
    j = series.iterations[ts]
    time = j.time / femto 
    E = j.meshes["E"]
    dz, dx = E.grid_spacing 
    Ez = E["z"]
    nz, nx = Ez.shape
    Ez_data = Ez.load_chunk()
    series.flush()
    gridx_w = dx*np.arange(0,nx)/micron        	
    gridy_w = dz*np.arange(0,nz)/micron        	
    lineoutx_w = Ez_data[int(nz*0.5),:]
    lineouty_w = Ez_data[:,int(nx*0.5)]
    print('warpx time [fs] = ', time)
    
    # smilei 
    Ey = s.Field(0,'Ey', units=["um", "V/m", "fs"])
    time = Ey.getTimes()[n]
    Ey_data = Ey.getData(ts)
    Ey_data = Ey_data[0]
    gridx_s = Ey.getAxis('x')
    gridy_s = Ey.getAxis('y')
    lineoutx_s = Ey_data[:,int(np.shape(Ey_data)[1]*0.5)]
    lineouty_s = Ey_data[int(np.shape(Ey_data)[0]*0.5),:]
    print('smilei time [fs] = ', time)
    
    # epoch
    fname = epoch_dir+'fields%04d.sdf' % n
    raw = sdf.read(fname)
    time = raw.Header["time"]/femto
    gridx_e = raw.__dict__["Grid_Grid_mid"].data[0]/micron
    gridy_e = raw.__dict__["Grid_Grid_mid"].data[1]/micron
    Ey = raw.__dict__["Electric_Field_Ey"].data
    lineoutx_e = Ey[:,int(np.shape(Ey)[1]*0.5)]
    lineouty_e = Ey[int(np.shape(Ey)[0]*0.5),:]
    print('epoch time [fs] = ', time)
        
    # plot images
    fig, ax = plt.subplots(ncols=2, nrows=1, figsize=(3000./my_dpi, 1000./my_dpi), dpi=my_dpi)

    ax[0].plot(gridx_w, lineoutx_w, lw=2, label='warpx')  
    ax[0].plot(gridx_s, lineoutx_s, lw=1, label='smilei')
    ax[0].plot(gridx_e, lineoutx_e, lw=1, label='epoch')

    ax[1].plot(gridy_w, lineouty_w, lw=1, label='warpx')  
    ax[1].plot(gridy_s, lineouty_s, lw=1, label='smilei')
    ax[1].plot(gridy_e, lineouty_e, lw=1, label='epoch')    
  
    ax[0].set_xlabel(r'x [$\mu$m]')
    ax[0].set_ylabel(r'E$_y$ [V/m]')
    ax[0].legend()
    ax[0].set_title('transverse field at mid y')



    #__________________________________________________________
    # save image
    image_file_name =plot_dir+'/lineouts_%04d.png' % ts
    plt.tight_layout()
    plt.savefig(image_file_name, dpi=my_dpi)# bbox_inches='tight', 
    plt.close()


