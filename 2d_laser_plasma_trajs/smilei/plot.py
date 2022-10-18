import numpy as np 
import happi
import matplotlib.pyplot as plt
from matplotlib.pyplot import * 
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import cm
from matplotlib.colors import LogNorm
from matplotlib.colors import Normalize

matplotlib.use('AGG')
my_dpi=300.
myunits = ['fs', 'um', 'MV/um', 'kT', 'J/um^2']


s=happi.Open('.')

micron = s.namelist.micron 
fs = s.namelist.fs
#_______________________________________________________________________________________________
# PLOT TRAJS ELE 

ele = s.TrackParticles(species='ele', sort=True, axes=['x', 'y', 'px', 'py', 'pz'])
Bz = s.Field(0, 'Bz', units=myunits)
Re = s.Field(0,'Rho_ele', units=myunits) 
x = Bz.getAxis('x')
y = Bz.getAxis('y') 

ts = ele.getAvailableTimesteps()

times = ele.getTimes()/fs

npart = 20
npart_tot=len(ele.getData().get('x')[0,:])

index = np.random.randint(0,npart_tot, npart)

ii = 1187500
i1 = int(ii*0.1)
i2 = int(ii*0.3)
 
index = np.argwhere(np.abs(ele.getData().get('y')[0,:]/micron) < 3)
index = np.random.randint(0,npart_tot, npart)
xtraj = ele.getData().get('x')[:,index]/micron
ytraj = ele.getData().get('y')[:,index]/micron


for n,i in enumerate(ts, start=0): 
    fig, ax = plt.subplots(1,3,figsize=(3000./my_dpi, 800./my_dpi), dpi=my_dpi)

    data = Re.getData(i)[0]
    im=ax[0].imshow(-np.transpose(data),alpha=1.0, vmin=0, vmax=1e-2,  extent=[np.min(x),np.max(x),np.min(y),np.max(y)], cmap='Greys')
    cbar = plt.colorbar(im, ax=ax[0])
    cbar.set_label(r'$n_e$ [$n_c$]') 


    data = Bz.getData(i)[0] 
    im=ax[1].imshow(np.transpose(data),extent=[np.min(x),np.max(x),np.min(y),np.max(y)], cmap='seismic', vmin=-50, vmax=50)
    cbar = plt.colorbar(im, ax=ax[1])
    cbar.set_label(r'B$_z$ [kT]')
    

    for p in range(npart):
        ax[2].plot(xtraj[0:n+1,p],ytraj[0:n+1,p], zorder=10, lw = 1)
        ax[2].scatter(xtraj[n,p],ytraj[n,p],zorder=11, s=1)

    for a in ax.reshape(-1):
        a.set_xlabel(r'x [$\mu$m]')
        a.set_ylabel(r'y [$\mu$m]')
        a.set_xlim(x.min(),x.max())
        a.axis('equal')

    plt.tight_layout() #rect=[0, 0.02, 1, 0.95])
    plt.savefig("FIG_%04d.png" % n)
    plt.close()

