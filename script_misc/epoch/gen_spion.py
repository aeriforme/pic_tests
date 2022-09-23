#!/galileo/home/userexternal/aforment/myenvnew/bin/python
import sdf_helper as sh 
from scipy.constants import * 
import sys
import numpy as np
import sdf 



lambda_SI = 0.8 * micron # wavelength
omega_SI = 2.0 * pi * c / lambda_SI
n_crit = 1.1 * 1e21 / centi**(3) * (1.e-6/lambda_SI)**(2) #m^-3
MeV_SI = 1e6*eV






path='/g100_scratch/userexternal/aforment/SCAN3D/'

list_dirs = ['DLT_t2_d1']

#list_dirs = ['DLT_t5_d1', 'DLT_t10_d1', 'DLT_t15_d1']
#list_dirs = ['NOFOAM', ]
#list_dirs= ['DLT_t20_d1', 'DLT_t25_d1']                                                                

for d in list_dirs:
    print(d) 
    data_dir = path + d + '/Data/'    
    filenames=np.genfromtxt(data_dir +'particles.visit',dtype='str')
    count = len(filenames)
    
    for t in range(count):

        fname = data_dir+'particles%04d.sdf' % t 
        raw=sdf.read(fname)
        time = raw.Header['time']/femto
        print('time [fs] = ', time)
        
        if hasattr(raw, 'Grid_Particles_ion_cont'):

            x=raw.__dict__["Grid_Particles_ion_cont"].data[0]/micron
            y=raw.__dict__["Grid_Particles_ion_cont"].data[1]/micron
            z=raw.__dict__["Grid_Particles_ion_cont"].data[2]/micron
            px=raw.__dict__["Particles_Px_ion_cont"].data/(m_p*c)
            py=raw.__dict__["Particles_Py_ion_cont"].data/(m_p*c)
            pz=raw.__dict__["Particles_Pz_ion_cont"].data/(m_p*c)
            w=raw.__dict__["Particles_Weight_ion_cont"].data


            E=m_p*c**2*(np.sqrt(1.+px**2+py**2+pz**2)-1.)/e*1e-6

            mybins = np.linspace(0,50,100) 
            H, b = np.histogram(E, bins = mybins, weights = w)
            head = "0) E [MeV] | 1) dN/dE [arb. units]"
            np.savetxt('SPION_E_dNdE_%s_%04d.txt' % (d, t), np.c_[b[:-1], H], header=head)
       
            







