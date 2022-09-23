import os 
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import micron,c,pi,centi,femto
import openpmd_api as io
from mpl_toolkits.axes_grid1 import make_axes_locatable
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
um = micron 
#__________________________________________________________
# simulation parameters
Lx = 70*micron 
Ly = 30*micron 
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
#__________________________________________________________
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
    #__________________________________________________________
    # plot images
    my_dpi = 300 
    fig, ax = plt.subplots(ncols=3, nrows=1, figsize=(3000./my_dpi, 1000./my_dpi), dpi=my_dpi, sharex=True)

    # warpx
    j = series.iterations[ts]
    time = j.time / femto 
    E = j.meshes["E"]
    Ez = E["z"]
    Ez_data = Ez.load_chunk()
    series.flush() 
    print('warpx shape = ', np.shape(Ez_data)) 	
    
    im = ax[0].imshow(Ez_data, cmap = 'seismic', extent = [0, Lx/um, 0, Ly/um])
    ax[0].set_title('warpx, t = %.2f [fs]' % time)
    divider = make_axes_locatable(ax[0])
    cax = divider.append_axes('right', size='5%', pad=0.05)
    fig.colorbar(im, cax=cax, orientation='vertical')

    # smilei 
    Ey = s.Field(0,'Ey', units=["um", "V/m", "fs"])
    time = Ey.getTimes()[n]
    Ey_data = Ey.getData(ts)
    Ey_data = Ey_data[0]
    print('smilei shape = ', np.shape(Ey_data)) 	

    im = ax[1].imshow(np.transpose(Ey_data), cmap ='seismic', extent = [0, Lx/um, 0, Ly/um])     
    ax[1].set_title('smilei, t = %.2f [fs]'% time)
    divider = make_axes_locatable(ax[1])
    cax = divider.append_axes('right', size='5%', pad=0.05)
    fig.colorbar(im, cax=cax, orientation='vertical')


    # epoch 
    fname = epoch_dir+'fields%04d.sdf' % n
    raw = sdf.read(fname)
    time = raw.Header["time"]/femto
    Ey = raw.__dict__["Electric_Field_Ey"].data
    print('epoch shape = ', np.shape(Ey)) 	

    im = ax[2].imshow(np.transpose(Ey), cmap ='seismic', extent = [0, Lx/um, 0, Ly/um])     
    ax[2].set_title('epoch, t = %.2f [fs]' % time)
    divider = make_axes_locatable(ax[2])
    cax = divider.append_axes('right', size='5%', pad=0.05)
    fig.colorbar(im, cax=cax, orientation='vertical')
  
  
    for a in ax.reshape(-1):
        a.set_xlabel(r'x [$\mu$m]')
    ax[0].set_ylabel(r'E$_y$ [V/m]')


    #__________________________________________________________
    # save image
    image_file_name =plot_dir+'/maps_%04d.png' % ts
    plt.tight_layout()
    plt.savefig(image_file_name,dpi=300)# bbox_inches='tight', 
    plt.close()


