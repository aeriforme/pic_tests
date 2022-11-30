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

#use('AGG')

#__________________________________________________________
# create plot directory
plot_dir = './plots'
if os.path.exists(plot_dir) is False:
    os.mkdir(plot_dir)

#__________________________________________________________
lambda_SI = 0.8*micron # wavelength
omega_SI = 2.0*pi*c / lambda_SI
fs = 1.e-15 * omega_SI;
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
series = io.Series(warpx_dir+"/openpmd_%T.h5",io.Access.read_only)
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
    fig, ax = plt.subplots(ncols=3, nrows=3, figsize=(4500./my_dpi, 3000./my_dpi), dpi=my_dpi, sharex=True)


    # smilei 
    ne_f = s.Field(0,'-Rho_ele_f')
    ne_s = s.Field(0,'-Rho_ele_s')
    ni_c = s.Field(0,'-Rho_ion_c')
    time = (ne_f.getTimes()[n])/fs
    data_s = ne_s.getData(ts)[0]
    data_f = ne_f.getData(ts)[0]
    data_c = -1*ni_c.getData(ts)[0]
    Bz = s.Field(0,'Bz', units=["um", "T", "fs"])
    data = Bz.getData(ts)[0] 
    bmax = data.max()
    bmin = data.min() 
    alphas = Normalize(0,bmax, clip=True)(np.abs(data))
    alphas = np.clip(alphas, 0.1, 1.) 
    colors = Normalize(bmin, bmax)(data)
    colors = cm.seismic(colors)
    colors[..., -1] = alphas
    
    ax[0,0].set_title('ele foam smilei, t = %.2f [fs]'% time)
    ax[0,1].set_title('ele subs smilei, t = %.2f [fs]'% time)
    ax[0,2].set_title('ion cont smilei, t = %.2f [fs]'% time)
    
    im=ax[0,0].imshow(np.transpose(data_f), vmin=0., vmax=0.2, extent=[0,Lx/um,0,Ly/um], cmap='Reds', aspect='equal')
    divider = make_axes_locatable(ax[0,0])
    cax = divider.append_axes('bottom', size='3%', pad=-0.5)
    fig.colorbar(im, cax=cax, orientation='horizontal', label=r'n$_e$ foam [n$_c$]')
    
    im=ax[0,0].imshow(np.transpose(colors,(1,0,2)), extent=[0,Lx/um,0,Ly/um], cmap='seismic', aspect='equal')
    divider = make_axes_locatable(ax[0,0])
    cax = divider.append_axes('bottom', size='3%', pad=1.1)
    fig.colorbar(im, cax=cax, orientation='horizontal', label=r'B$_z$ [T]')
    
    im=ax[0,1].imshow(np.transpose(data_s), vmin=0., vmax=50, extent=[0,Lx/um,0,Ly/um], cmap='binary', aspect='equal')
    divider = make_axes_locatable(ax[0,1])
    cax = divider.append_axes('bottom', size='3%', pad=-0.5)
    fig.colorbar(im, cax=cax, orientation='horizontal', label=r'n$_e$ subs [n$_c$]')
    
    im=ax[0,1].imshow(np.transpose(colors,(1,0,2)), extent=[0,Lx/um,0,Ly/um], cmap='seismic', aspect='equal')
    divider = make_axes_locatable(ax[0,1])
    cax = divider.append_axes('bottom', size='3%', pad=1.1)
    fig.colorbar(im, cax=cax, orientation='horizontal', label=r'B$_z$ [T]')
    
    im=ax[0,2].imshow(np.transpose(data_c), vmin=0., vmax=50, extent=[0,Lx/um,0,Ly/um], cmap='Blues', aspect='equal')
    divider = make_axes_locatable(ax[0,2])
    cax = divider.append_axes('bottom', size='3%', pad=-0.5)
    fig.colorbar(im, cax=cax, orientation='horizontal', label=r'n$_i$ cont [n$_c$]')
    
    im=ax[0,2].imshow(np.transpose(colors,(1,0,2)), extent=[0,Lx/um,0,Ly/um], cmap='seismic', aspect='equal')
    divider = make_axes_locatable(ax[0,2])
    cax = divider.append_axes('bottom', size='3%', pad=1.1)
    fig.colorbar(im, cax=cax, orientation='horizontal', label=r'B$_z$ [T]')    
    



    # epoch 
    fname = epoch_dir+'fields%04d.sdf' % n
    raw = sdf.read(fname)
    time = raw.Header["time"]/femto
    ne_f = raw.__dict__["Derived_Number_Density_ele_foam"].data
    ne_s = raw.__dict__["Derived_Number_Density_ele_subs"].data 
    ni_c = raw.__dict__["Derived_Number_Density_ion_cont"].data 
    data_f = ne_f / n_crit
    data_s = ne_s / n_crit
    data_c = ni_c / n_crit
    Bz = raw.__dict__["Magnetic_Field_Bz"].data
    data = Bz
    bmax = data.max()
    bmin = data.min() 
    alphas = Normalize(0,bmax, clip=True)(np.abs(data))
    alphas = np.clip(alphas, 0.1, 1.) 
    colors = Normalize(bmin, bmax)(data)
    colors = cm.seismic(colors)
    colors[..., -1] = alphas
    
    ax[1,0].set_title('ele foam epoch, t = %.2f [fs]' % time)
    ax[1,1].set_title('ele subs epoch, t = %.2f [fs]' % time)
    ax[1,2].set_title('ion cont epoch, t = %.2f [fs]' % time)
  
    im=ax[1,0].imshow(np.transpose(data_f), vmin=0., vmax=0.2, extent=[0,Lx/um,0,Ly/um], cmap='Reds', aspect='equal')
    divider = make_axes_locatable(ax[1,0])
    cax = divider.append_axes('bottom', size='3%', pad=-0.5)
    fig.colorbar(im, cax=cax, orientation='horizontal', label=r'n$_e$ foam [n$_c$]')
    
    im=ax[1,0].imshow(np.transpose(colors,(1,0,2)), extent=[0,Lx/um,0,Ly/um], cmap='seismic', aspect='equal')
    divider = make_axes_locatable(ax[1,0])
    cax = divider.append_axes('bottom', size='3%', pad=1.1)
    fig.colorbar(im, cax=cax, orientation='horizontal', label=r'B$_z$ [T]')
    
    im=ax[1,1].imshow(np.transpose(data_s), vmin=0., vmax=50, extent=[0,Lx/um,0,Ly/um], cmap='binary', aspect='equal')
    divider = make_axes_locatable(ax[1,1])
    cax = divider.append_axes('bottom', size='3%', pad=-0.5)
    fig.colorbar(im, cax=cax, orientation='horizontal', label=r'n$_e$ subs [n$_c$]')
    
    im=ax[1,1].imshow(np.transpose(colors,(1,0,2)), extent=[0,Lx/um,0,Ly/um], cmap='seismic', aspect='equal')
    divider = make_axes_locatable(ax[1,1])
    cax = divider.append_axes('bottom', size='3%', pad=1.1)
    fig.colorbar(im, cax=cax, orientation='horizontal', label=r'B$_z$ [T]')
    
    im=ax[1,2].imshow(np.transpose(data_c), vmin=0., vmax=50, extent=[0,Lx/um,0,Ly/um], cmap='Blues', aspect='equal')
    divider = make_axes_locatable(ax[1,2])
    cax = divider.append_axes('bottom', size='3%', pad=-0.5)
    fig.colorbar(im, cax=cax, orientation='horizontal', label=r'n$_i$ cont [n$_c$]')
    
    im=ax[1,2].imshow(np.transpose(colors,(1,0,2)), extent=[0,Lx/um,0,Ly/um], cmap='seismic', aspect='equal')
    divider = make_axes_locatable(ax[1,2])
    cax = divider.append_axes('bottom', size='3%', pad=1.1)
    fig.colorbar(im, cax=cax, orientation='horizontal', label=r'B$_z$ [T]')       
  


    # warpx
    j = series.iterations[ts]
    time = j.time / femto       
    rho_ele_f = j.meshes["rho_elef"][io.Mesh_Record_Component.SCALAR]
    rho_ele_s = j.meshes["rho_eles"][io.Mesh_Record_Component.SCALAR]
    rho_ion_c = j.meshes["rho_ionc"][io.Mesh_Record_Component.SCALAR]
    rho_ele_data_f = rho_ele_f.load_chunk()
    series.flush() 
    data_f = rho_ele_data_f / (-e * n_crit) 
    rho_ele_data_s = rho_ele_s.load_chunk()
    series.flush() 
    data_s = rho_ele_data_s / (-e * n_crit) 
    rho_ion_data_c = rho_ion_c.load_chunk()
    series.flush() 
    data_c = rho_ion_data_c / (e * n_crit) 
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

    ax[2,0].set_title('ele foam warpx, t = %.2f [fs]' % time)
    ax[2,1].set_title('ele subs warpx, t = %.2f [fs]' % time)
    ax[2,2].set_title('ion cont warpx, t = %.2f [fs]' % time)

    im=ax[2,0].imshow(data_f, vmin=0., vmax=0.2, extent=[0,Lx/um,0,Ly/um], cmap='Reds', aspect='equal')
    divider = make_axes_locatable(ax[2,0])
    cax = divider.append_axes('bottom', size='3%', pad=-0.5)
    fig.colorbar(im, cax=cax, orientation='horizontal', label=r'n$_e$ foam [n$_c$]')
    
    im=ax[2,0].imshow(colors, extent=[0,Lx/um,0,Ly/um], cmap='seismic', aspect='equal')
    divider = make_axes_locatable(ax[2,0])
    cax = divider.append_axes('bottom', size='3%', pad=1.1)
    fig.colorbar(im, cax=cax, orientation='horizontal', label=r'B$_z$ [T]')
    
    im=ax[2,1].imshow(data_s, vmin=0., vmax=50, extent=[0,Lx/um,0,Ly/um], cmap='binary', aspect='equal')
    divider = make_axes_locatable(ax[2,1])
    cax = divider.append_axes('bottom', size='3%', pad=-0.5)
    fig.colorbar(im, cax=cax, orientation='horizontal', label=r'n$_e$ subs [n$_c$]')
    
    im=ax[2,1].imshow(colors, extent=[0,Lx/um,0,Ly/um], cmap='seismic', aspect='equal')
    divider = make_axes_locatable(ax[2,1])
    cax = divider.append_axes('bottom', size='3%', pad=1.1)
    fig.colorbar(im, cax=cax, orientation='horizontal', label=r'B$_z$ [T]')
    
    im=ax[2,2].imshow(data_c, vmin=0., vmax=50, extent=[0,Lx/um,0,Ly/um], cmap='Blues', aspect='equal')
    divider = make_axes_locatable(ax[2,2])
    cax = divider.append_axes('bottom', size='3%', pad=-0.5)
    fig.colorbar(im, cax=cax, orientation='horizontal', label=r'n$_i$ cont [n$_c$]')
    
    im=ax[2,2].imshow(colors, extent=[0,Lx/um,0,Ly/um], cmap='seismic', aspect='equal')
    divider = make_axes_locatable(ax[2,2])
    cax = divider.append_axes('bottom', size='3%', pad=1.1)
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


