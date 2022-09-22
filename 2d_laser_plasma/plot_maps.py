import os 
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import micron,c,pi,centi,femto,e
import openpmd_api as io
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import LogNorm, Normalize
import happi 
import sdf 
from matplotlib import use, cm

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
    ax[0].set_title('warpx, t = %.2f [fs]' % time)
        
    rho_ele = j.meshes["rho_ele"][io.Mesh_Record_Component.SCALAR]
    rho_ele_data = rho_ele.load_chunk()
    series.flush() 
    data = rho_ele_data / (-e * n_crit) 
    im=ax[0].imshow(data, vmin=0., vmax=0.02, extent=[0,Lx/um,0,Ly/um], cmap='binary', aspect='equal')
    divider = make_axes_locatable(ax[0])
    cax = divider.append_axes('bottom', size='5%', pad=0.02)
    fig.colorbar(im, cax=cax, orientation='horizontal',label=r'n$_e$ [n$_c$]')

    B = j.meshes["B"]
    Bz = B["y"]
    data = Bz.load_chunk()
    series.flush() 
    bmax = data.max()
    bmin = data.min() 
    alphas = Normalize(0,bmax, clip=True)(np.abs(data))
    alphas = np.clip(alphas, 0.1, 1.) 
    colors = Normalize(bmin, bmax)(data)
    colors = cm.seismic(colors)
    colors[..., -1] = alphas

    im=ax[0].imshow(colors, extent=[0,Lx/um,0,Ly/um], cmap='seismic', aspect='equal')
    ax[0].set_title('smilei, t = %.2f [fs]'% time)
    divider = make_axes_locatable(ax[0])
    cax = divider.append_axes('bottom', size='5%', pad=1.2)
    fig.colorbar(im, cax=cax, orientation='horizontal', label=r'B$_z$ [T]')


    # smilei 
    ne = s.Field(0,'-Rho_ele')
    time = ne.getTimes()[n]
    ax[1].set_title('smilei, t = %.2f [fs]'% time)
    
    data = ne.getData(ts)[0]
    im=ax[1].imshow(np.transpose(data), vmin=0., vmax=0.02, extent=[0,Lx/um,0,Ly/um], cmap='binary', aspect='equal')
    divider = make_axes_locatable(ax[1])
    cax = divider.append_axes('bottom', size='5%', pad=0.02)
    fig.colorbar(im, cax=cax, orientation='horizontal', label=r'n$_e$ [n$_c$]')


    Bz = s.Field(0,'Bz', units=["um", "T", "fs"])
    data = Bz.getData(ts)[0] 
    bmax = data.max()
    bmin = data.min() 
    alphas = Normalize(0,bmax, clip=True)(np.abs(data))
    alphas = np.clip(alphas, 0.1, 1.) 
    colors = Normalize(bmin, bmax)(data)
    colors = cm.seismic(colors)
    colors[..., -1] = alphas
    im=ax[1].imshow(np.transpose(colors,(1,0,2)), extent=[0,Lx/um,0,Ly/um], cmap='seismic', aspect='equal')
    divider = make_axes_locatable(ax[1])
    cax = divider.append_axes('bottom', size='5%', pad=1.2)
    fig.colorbar(im, cax=cax, orientation='horizontal', label=r'B$_z$ [T]')


    # epoch 
    fname = epoch_dir+'fields%04d.sdf' % n
    raw = sdf.read(fname)
    time = raw.Header["time"]/femto
    ax[2].set_title('epoch, t = %.2f [fs]' % time)

    ne = raw.__dict__["Derived_Number_Density_ele"].data 
    data = ne / n_crit 
    
    im=ax[2].imshow(np.transpose(data), vmin=0., vmax=0.02, extent=[0,Lx/um,0,Ly/um], cmap='binary', aspect='equal')
    divider = make_axes_locatable(ax[2])
    cax = divider.append_axes('bottom', size='5%', pad=0.02)
    fig.colorbar(im, cax=cax, orientation='horizontal', label=r'n$_e$ [n$_c$]')

    Bz = raw.__dict__["Magnetic_Field_Bz"].data
    data = Bz
    bmax = data.max()
    bmin = data.min() 
    alphas = Normalize(0,bmax, clip=True)(np.abs(data))
    alphas = np.clip(alphas, 0.1, 1.) 
    colors = Normalize(bmin, bmax)(data)
    colors = cm.seismic(colors)
    colors[..., -1] = alphas
    im=ax[2].imshow(np.transpose(colors,(1,0,2)), extent=[0,Lx/um,0,Ly/um], cmap='seismic', aspect='equal')
    divider = make_axes_locatable(ax[2])
    cax = divider.append_axes('bottom', size='5%', pad=1.2)
    fig.colorbar(im, cax=cax, orientation='horizontal', label=r'B$_z$ [T]')

    
  
  
    for a in ax.reshape(-1):
        a.set_xlabel(r'x [$\mu$m]')
        a.set_ylabel(r'y [$\mu$m]')


    #__________________________________________________________
    # save image
    image_file_name =plot_dir+'/maps_%04d.png' % ts
    plt.tight_layout()
    plt.savefig(image_file_name,dpi=300)# bbox_inches='tight', 
    plt.close()


