import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import micron,c,pi,centi,femto,elementary_charge as q_e
import openpmd_api as io
import happi 
import sdf 
from matplotlib import use

#use('AGG')

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
Lx = 70.*micron 
resx = 20. # points per micron 
nx = Lx * resx / micron 
dx = Lx/nx

cfl = 0.98
dt = cfl * dx / (c * np.sqrt(2.))

Tsim = 1.5*Lx / c 
nsteps = int(Tsim/dt) 


every_fs = np.floor(femto/dt)
out_freq = int(10*every_fs)

steps = np.arange(0,nsteps,out_freq)

# warpx
warpx_dir='./warpx/diags/Fields'
# series to know where the data is
series = io.Series(warpx_dir+"/openpmd_%T.bp",io.Access.read_only)
#If you only had warpx, n and ts could have been taken with the following:
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
    ################
    # warpx
    j = series.iterations[ts]
    time_w = j.time / femto 
    #extracting E
    E = j.meshes["E"]
    dz, dx = E.grid_spacing 
    Ez = E["z"]
    nz, nx = Ez.shape
    Ez_data = Ez.load_chunk()
    
    #extracting B
    B = j.meshes["B"]
    By = B["y"]
    By_data = By.load_chunk()
    
    #extracting Rho ele
    Rhoe = j.meshes["rho_ele"][io.Mesh_Record_Component.SCALAR]
    Rhoe_data = Rhoe.load_chunk()
    
    #actual loading of data
    series.flush()    
    
    #making the grid 
    gridx_w = dx*np.arange(0,nx)/micron        	
    gridy_w = dz*np.arange(0,nz)/micron        	
    
    lineout_Ez_w = Ez_data[int(nz*0.5),:]
    lineout_By_w = By_data[int(nz*0.5),:]
    lineout_Rhoe_w = Rhoe_data[int(nz*0.5),:]
    lineout_ne_w = lineout_Rhoe_w/(-q_e*n_crit)
    
    print('warpx shape = ', np.shape(Ez_data))
    print('warpx max = ', np.max(Ez_data)/1e12)
    
    ###############
    # smilei 
    Ey = s.Field(0,'Ey', units=["um", "V/m", "fs"], timesteps=ts)
    Bz = s.Field(0,'Bz', units=["um", "T", "fs"], timesteps=ts)
    ne = s.Field(0,'Rho_ele', units=["um", "fs"], timesteps=ts)
    
    time_s = Ey.getTimes()
    
    Ey_data = Ey.getData()
    Bz_data = Bz.getData()
    ne_data = ne.getData()   
    # actual data in the first element of object xx.getData(ts)
    Ey_data = Ey_data[0]
    Bz_data = Bz_data[0]
    ne_data = ne_data[0]
    
    gridx_s = Ey.getAxis('x')
    gridy_s = Ey.getAxis('y')
    
    lineout_Ey_s = Ey_data[:,int(np.shape(Ey_data)[1]*0.5)]
    lineout_Bz_s = Bz_data[:,int(np.shape(Bz_data)[1]*0.5)]
    lineout_ne_s = ne_data[:,int(np.shape(ne_data)[1]*0.5)]


    print('smilei shape = ', np.shape(Ey_data))
    print('smilei max = ', np.max(Ey_data)/1e12)

    ################
    # epoch
    fname = epoch_dir+'fields%04d.sdf' % n  # fname = filenames[n]
    raw = sdf.read(fname)
    time_e = raw.Header["time"]/femto
    gridx_e = raw.__dict__["Grid_Grid_mid"].data[0]/micron
    gridy_e = raw.__dict__["Grid_Grid_mid"].data[1]/micron
    Ey = raw.__dict__["Electric_Field_Ey"].data
    Bz = raw.__dict__["Magnetic_Field_Bz"].data
    ne = (raw.__dict__["Derived_Number_Density_ele"].data)/n_crit
    
    lineout_Ey_e = Ey[:,int(np.shape(Ey)[1]*0.5)]
    lineout_Bz_e = Bz[:,int(np.shape(Bz)[1]*0.5)]
    lineout_ne_e = ne[:,int(np.shape(ne)[1]*0.5)]

    ##################
    # plot images
    fig, ax = plt.subplots(ncols=3, nrows=1, figsize=(3000./my_dpi, 1000./my_dpi), dpi=my_dpi)

    ax[0].plot(gridx_w, lineout_Ez_w, lw=1, label='warpx')  
    ax[0].plot(gridx_s, lineout_Ey_s, lw=1, label='smilei')
    ax[0].plot(gridx_e, lineout_Ey_e, lw=1, label='epoch')

    ax[1].plot(gridx_w, lineout_By_w, lw=1, label='warpx')  
    ax[1].plot(gridx_s, lineout_Bz_s, lw=1, label='smilei')
    ax[1].plot(gridx_e, lineout_Bz_e, lw=1, label='epoch') 
    
    ax[2].plot(gridx_w, lineout_ne_w, lw=1, label='warpx')  
    ax[2].plot(gridx_s, lineout_ne_s, lw=1, label='smilei')
    ax[2].plot(gridx_e, lineout_ne_e, lw=1, label='epoch')     
  
    ax[0].set_xlabel(r'x [$\mu$m]')
    ax[0].set_ylabel(r'E$_y$ [V/m]')
    ax[0].legend()
    ax[0].set_title('transverse E field at mid y')
    
    ax[1].set_xlabel(r'x [$\mu$m]')
    ax[1].set_ylabel(r'B$_z$ [T]')
    ax[1].legend()
    ax[1].set_title('transverse B field at mid y')

    ax[2].set_xlabel(r'x [$\mu$m]')
    ax[2].set_ylabel(r'n$_e$ [n$_{crit}$]')
    ax[2].legend()
    ax[2].set_title('electron density at mid y')


    #__________________________________________________________
    # save image
    image_file_name =plot_dir+'/lineouts_%04d.png' % ts
    plt.tight_layout()
    plt.savefig(image_file_name, dpi=my_dpi)# bbox_inches='tight', 
    plt.close()


