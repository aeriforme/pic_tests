#!/galileo/home/userexternal/aforment/myenvnew/bin/python
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
norm = LogNorm(vmin=0.1,vmax=1e5,clip=False)
norm = Normalize(vmin=0.,vmax=1e5,clip=False)
norm = LogNorm(vmin=1e4,vmax=1e7) #,clip=False)


path='/g100_scratch/userexternal/aforment/SCAN3D/'
list_dirs = ['DLT_t2_d1', 'DLT_t5_d1', 'DLT_t10_d1', 'DLT_t15_d1', 'DLT_t20_d1', 'DLT_t25_d1']
#list_dirs = 


bin_phi =np.linspace(-np.pi,np.pi,200)
bin_theta = np.linspace(-np.pi,np.pi,200)

for d in list_dirs:
    print(d) 
    data_dir = path + d + '/Data/'    
    filenames=np.genfromtxt(data_dir +'particles.visit',dtype='str')
    count = len(filenames)
    
    for t in range(count-1, count):
        fname = data_dir+'particles%04d.sdf' % t 
        raw=sdf.read(fname)
        #time = raw.Header['time']/femto
        #print('time [fs] = ', time)
        
        if hasattr(raw, 'Grid_Particles_pho'):

            fig, ax=plt.subplots(1,4,dpi=300, figsize=(12.5,4), gridspec_kw={"width_ratios":[1,1,1, 0.05]})

            x=raw.__dict__["Grid_Particles_pho"].data[0]/micron
            y=raw.__dict__["Grid_Particles_pho"].data[1]/micron
            z=raw.__dict__["Grid_Particles_pho"].data[2]/micron
            px=raw.__dict__["Particles_Px_pho"].data
            py=raw.__dict__["Particles_Py_pho"].data
            pz=raw.__dict__["Particles_Pz_pho"].data

            pp = np.sqrt(px**2+py**2+pz**2)
            E=c*pp/e*1e-6

            phi=np.arctan2(py,px) #*180./np.pi
            theta=np.arccos(pz/pp) #*180./np.pi

            w = raw.__dict__["Particles_Weight_pho"].data

            # 1 MeV < E < 5 MeV
            H, bx, by = np.histogram2d(phi[(E>1.)&(E<5.)],theta[(E>1.)&(E<5.)], bins=(bin_phi, bin_theta), weights=w[(E>1.)&(E<5.)])            
            print(H.min(), H.max())

            ax[0].imshow(np.transpose(H),extent=[bin_phi[0], bin_phi[-1], bin_theta[0], bin_theta[-1]], norm=norm) # , norm=LogNorm(vmin=H.min(), vmax=H.max()))
            ax[0].set_aspect('equal')

            # 5 MeV < E < 10 MeV
            H, bx, by = np.histogram2d(phi[(E>5.)&(E<10.)],theta[(E>5.)&(E<10.)], bins=(bin_phi, bin_theta), weights=w[(E>5.)&(E<10.)])   	      	
            ax[1].imshow(np.transpose(H),extent=[bin_phi[0], bin_phi[-1], bin_theta[0], bin_theta[-1]], norm=norm) # norm=LogNorm(vmin=H.min(), vmax=H.max()) )
            ax[1].set_aspect('equal')
            print(H.min(), H.max())
            
            #  E > 10 MeV                                                                                                                                                        
            H, bx, by = np.histogram2d(phi[E>10.],theta[E>10.], bins=(bin_phi, bin_theta), weights=w[E>10.])
            im = ax[2].imshow(np.transpose(H),extent=[bin_phi[0], bin_phi[-1], bin_theta[0], bin_theta[-1]], norm=norm) #LogNorm(vmin=H.min(), vmax=H.max()))
            ax[2].set_aspect('equal')
            print(H.min(), H.max())


            for a in ax.reshape(-1):
                a.set_xlabel(r'$\phi$')
                a.set_ylabel(r'$\theta$')

            fig.colorbar(im, cax=ax[3], label = 'xxx [arb. units]')
                
            plt.tight_layout()
            plt.savefig("rcf_pho_%s_%04d.png" % (d,t))
            plt.close('all')
            
            #np.savetxt('PHO2_X_Y_Z_PX_PY_PZ_E_%04d_%s.txt' % (t,d), np.c_[x,y,z,px,py,pz,E]) #,fmt='%.6e')







