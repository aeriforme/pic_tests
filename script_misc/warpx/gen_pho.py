import openpmd_api as io
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import micron, centi, pi, c, e, m_e, eV 
import vtk
from vtk.util import numpy_support
import re

MeV = 1e6*eV
lambda_SI = 0.8 * micron # wavelength
omega_SI = 2.0 * pi * c / lambda_SI
n_crit = 1.1 * 1e21 / centi**(3) * (1.e-6/lambda_SI)**(2) #m^-3

rootm = '/m100_scratch/userexternal/mgalbiat/3DSCANNICS/'
roota = '/m100_scratch/userexternal/aforment/3DSCANSYNCHR/'

dirs = (rootm+'DLT_a020_l0_d1', rootm+'DLT_a020_l2_d1', rootm+'DLT_a020_l5_d1', roota+'DLT_a020_l10_d1', roota+'DLT_a020_l15_d1', roota+'DLT_a020_l20_d1', roota+'DLT_a020_l25_d1') 

#names = ('l0', 'l2', 'l5', 'l10', 'l15', 'l20', 'l25') 

names = []
for k,d in enumerate(dirs):

    names.append( re.search(r'DLT_a020_l(.*?)_d1', d).group(1) ) 
    print(k,d,names[k]) 
    
    path =  d + '/diags/Fields_Particles/'

    series = io.Series(path+"openpmd_%T.bp", io.Access.read_only)
    #print("The Series contains {0} iterations:".format(len(series.iterations)))


    for n,i in enumerate(series.iterations,start=0):
        print("Iteration number \t {0}".format(i))
        j = series.iterations[i]

        #print("Current iteration contains {0} meshes:".format(len(j.meshes)))
        #for m in j.meshes:
        #    print("\t {0}".format(m))
        #    print("")
        #    print("Current iteration contains {0} particle species:".format(len(j.particles)))
        #     for ps in j.particles:
        #          print("\t {0}".format(ps))
        #         print("With records:")
        #         for r in j.particles[ps]:
        #             print("\t {0}".format(r))

        X, Y, Z, PX, PY, PZ, E, W = [], [], [], [], [], [], [], []

        if names[k]=='0': 
            species = ("pho_subs", "pho_cont")
        else:
            species = ("pho_foam", "pho_subs", "pho_cont")
            
        for spec in species:
        
            pho = j.particles[spec]
            x = pho["position"]["x"]
            y = pho["position"]["y"]
            z = pho["position"]["z"]
            px = pho["momentum"]["x"]
            py = pho["momentum"]["y"]
            pz = pho["momentum"]["z"]
            w = pho["weighting"][io.Mesh_Record_Component.SCALAR]

            x_data = x.load_chunk()
            y_data = y.load_chunk()
            z_data = z.load_chunk()
            px_data = px.load_chunk()
            py_data = py.load_chunk()
            pz_data = pz.load_chunk()   
            w_data = w.load_chunk()
            
            series.flush()

            en = np.sqrt(px_data**2+py_data**2+pz_data**2)*c

            X=np.append(X,x_data)
            Y=np.append(Y,y_data)
            Z=np.append(Z,z_data)
            PX=np.append(PX,px_data)
            PY=np.append(PY,py_data)
            PZ=np.append(PZ,pz_data)
            E=np.append(E,en)
            W=np.append(W,w_data)
            
        head='#0) x[um] | 1) y[um] | 2) z[um] | 3) px [m_e*c] | 4) py [m_e*c] | 5) pz [m_e*c] | 6) E [MeV] | 7) w [-] '
        np.savetxt('PHO_X_Y_Z_PX_PY_PZ_E_W_l%s_%04d.txt' % (names[k], i), np.c_[X/micron,Y/micron,Z/micron,PX/(m_e*c),PY/(m_e*c),PZ/(m_e*c),E/MeV, W],fmt='%.6e', header=head)

    del series 
