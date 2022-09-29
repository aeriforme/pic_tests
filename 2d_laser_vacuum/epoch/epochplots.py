import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import micron,c,pi,centi,femto
import sdf 
from matplotlib import use
from mpl_toolkits.axes_grid1 import make_axes_locatable

use('AGG') # per evitare che si aprano mille finestre

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
my_dpi = 300 


#__________________________________________________________
# simulation parameters
Lx = 70*micron 
resx = 25. # points per micron 
nx = Lx * resx / micron 
dx = Lx/nx

Ly = 30*micron

cfl = 0.98
dt = cfl * dx / (c * np.sqrt(2.))

Tsim = 1.5*Lx / c 
nsteps = int(Tsim/dt) 


every_fs = np.floor(femto/dt)
out_freq = int(10*every_fs)

steps = out_freq * np.arange(nsteps)
times = steps * dt/femto

steps = np.arange(0,nsteps,out_freq)

# epoch
epoch_dir = './Data/'
filenames=np.genfromtxt(epoch_dir+'fields.visit',dtype='str') 
#for t in range(len(filenames)): 

for n,ts in enumerate(steps):
   
    # epoch
    fname = epoch_dir+'fields%04d.sdf' % n
    raw = sdf.read(fname)
    time = raw.Header["time"]/femto
    gridx = raw.__dict__["Grid_Grid_mid"].data[0]/micron
    gridy = raw.__dict__["Grid_Grid_mid"].data[1]/micron
    Ey = raw.__dict__["Electric_Field_Ey"].data
    lineoutx = Ey[:,int(np.shape(Ey)[1]*0.5)]
    lineouty = Ey[int(np.shape(Ey)[0]*0.5),:]
    print('epoch time [fs] = ', time)
        
    # plot images
    fig, ax = plt.subplots(ncols=3, nrows=1, figsize=(4000./my_dpi, 1200./my_dpi), dpi=my_dpi)

    ax[0].plot(gridx, lineoutx, lw=1)
    ax[0].set_xlabel(r'x [$\mu$m]')
    ax[0].set_ylabel(r'E$_y$ [V/m]')
    ax[0].set_title('transverse field at mid y')

    ax[1].plot(gridy, lineouty, lw=1)    
    ax[1].set_xlabel(r'y [$\mu$m]')
    ax[1].set_ylabel(r'E$_y$ [V/m]')
    ax[1].set_title('transverse field at mid x')

    im = ax[2].imshow(np.transpose(Ey), cmap ='seismic', extent = [0, Lx/um, 0, Ly/um])     
    ax[2].set_xlabel(r'x [$\mu$m]')
    ax[2].set_ylabel(r'y [$\mu$m]')
    ax[2].set_title(r'E$_y$ [V/m] map')
    divider = make_axes_locatable(ax[2])
    cax = divider.append_axes('right', size='5%', pad=0.05)
    fig.colorbar(im, cax=cax, orientation='vertical')
  
    #__________________________________________________________
    # save image
    image_file_name =plot_dir+'/fig_%04d.png' % ts
    plt.tight_layout()
    fig.suptitle('epoch, t = %.2f [fs]' % time)

    plt.savefig(image_file_name, dpi=my_dpi)# bbox_inches='tight', 
    plt.close()


