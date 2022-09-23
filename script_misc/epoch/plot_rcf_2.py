import sdf_helper as sh 
from scipy.constants import * 
import sys
import numpy as np
import sdf 
import matplotlib.pyplot as plt 


lambda_SI = 0.8 * micron # wavelength
omega_SI = 2.0 * pi * c / lambda_SI
n_crit = 1.1 * 1e21 / centi**(3) * (1.e-6/lambda_SI)**(2) #m^-3
MeV_SI = 1e6*eV


from matplotlib.colors import ListedColormap, BoundaryNorm, Normalize, LogNorm
#norm = LogNorm(vmin=0.1,vmax=1e5,clip=False)
#norm = Normalize(vmin=0.,vmax=1e5,clip=False)
#norm = LogNorm(vmin=1e2,vmax=1e4,clip=False)

norm = LogNorm(vmin=1e-2,vmax=1,clip=False)


Emin = 1
Emax = 100


path='/g100_scratch/userexternal/aforment/SCAN3D/'
list_dirs = ['NOFOAM', 'DLT_t2_d1', 'DLT_t5_d1', 'DLT_t10_d1', 'DLT_t15_d1', 'DLT_t20_d1', 'DLT_t25_d1']
lab = [r't = 0 $\mu$m', r't = 2 $\mu$m', r't = 5 $\mu$m',  r't = 10 $\mu$m', r't = 15 $\mu$m',r't = 20 $\mu$m',r't = 25 $\mu$m']


bin_phi =np.linspace(-90.,90.,200)
bin_theta = np.linspace(-90.,90.,200)

fig, ax=plt.subplots(1,8,dpi=300, figsize=(24,4.2), gridspec_kw={"width_ratios":[1,1,1,1,1,1,1, 0.05]})

for i,d in enumerate(list_dirs, start=0):
    print(i, d) 
    data_dir = path + d + '/Data/'    
    filenames=np.genfromtxt(data_dir +'particles.visit',dtype='str')
    count = len(filenames)
    t = count-1
        
    fname = data_dir+'particles%04d.sdf' % t 
    raw=sdf.read(fname)
            
    if hasattr(raw, 'Grid_Particles_pho'):
        
        x=raw.__dict__["Grid_Particles_pho"].data[0]/micron
        y=raw.__dict__["Grid_Particles_pho"].data[1]/micron
        z=raw.__dict__["Grid_Particles_pho"].data[2]/micron
        px=raw.__dict__["Particles_Px_pho"].data
        py=raw.__dict__["Particles_Py_pho"].data
        pz=raw.__dict__["Particles_Pz_pho"].data

        E=c*np.sqrt(px**2+py**2+pz**2)/e*1e-6

        phi=np.arctan2(py,px)*180./np.pi

        rr = np.sqrt(px**2+py**2)
        theta=np.arctan2(pz,rr)*180./np.pi

        w = raw.__dict__["Particles_Weight_pho"].data
            
        index = (E>Emin) & (E<Emax)
        H, bx, by = np.histogram2d(phi[index],theta[index], bins=(bin_phi, bin_theta), weights=w[index], density=True)
        im = ax[i].imshow(np.transpose(H),extent=[bin_phi[0], bin_phi[-1], bin_theta[0], bin_theta[-1]], cmap='turbo', interpolation='nearest') # norm = norm  
        ax[i].set_aspect('equal')
        print(H.min(), H.max())
        ax[i].set_title(lab[i])


for a in ax.reshape(-1)[:-1]:
    a.set_xlabel(r'$\phi$')
    
for a in ax.reshape(-1)[1:-1]:    
    a.yaxis.set_ticklabels([])

ax[0].set_ylabel(r'$\theta$')
fig.colorbar(im, cax=ax[-1], label = r'd$^2$N/d$\phi$d$\theta$ [arb. units]')                

plt.subplots_adjust( left = 0.0, right = 0.96 , hspace=0.0, wspace=0.1)
#plt.subplots_adjust(wspace=0., hspace=0.)
plt.tight_layout()
plt.savefig("rcf_pho_Emin_%g_Emax_%g_%04d.png" % (Emin, Emax, t) )
plt.close('all')
            
#np.savetxt('PHO2_X_Y_Z_PX_PY_PZ_E_%04d_%s.txt' % (t,d), np.c_[x,y,z,px,py,pz,E]) #,fmt='%.6e')







