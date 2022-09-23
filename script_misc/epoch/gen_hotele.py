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

#list_dirs = ['NOFOAM', 'DLT_t2_d1', 'DLT_t5_d1', 'DLT_t10_d1', 'DLT_t15_d1', 'DLT_t20_d1', 'DLT_t25_d1']

list_dirs = ['DLT_t15_d1', ]

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
        
        if hasattr(raw, 'Grid_Particles_ele_foam'):

            x=raw.__dict__["Grid_Particles_ele_foam"].data[0]/micron
            y=raw.__dict__["Grid_Particles_ele_foam"].data[1]/micron
            z=raw.__dict__["Grid_Particles_ele_foam"].data[2]/micron
            px=raw.__dict__["Particles_Px_ele_foam"].data/(m_e*c)
            py=raw.__dict__["Particles_Py_ele_foam"].data/(m_e*c)
            pz=raw.__dict__["Particles_Pz_ele_foam"].data/(m_e*c)

            E=m_e*c**2*(np.sqrt(1.+px**2+py**2+pz**2)-1.)/e*1e-6

            Emin=10 #MeV
            
            print('number of macro-photons = ', len(E))
            #index = np.random.randint(0, high=len(E), size=int(0.5*len(E)))
            #np.savetxt('PHO_X_Y_Z_PX_PY_PZ_E_%04d_%s.txt' % (t,d), np.c_[x[index],y[index],z[index],px[index],py[index],pz[index],E[index]],fmt='%.6e')

            
            index = (E>=Emin) 
            np.savetxt('ELEFOAM_X_Y_Z_PX_PY_PZ_E_%04d_%s.txt' % (t,d), np.c_[x[index],y[index],z[index],px[index],py[index],pz[index],E[index]],fmt='%.6e')







