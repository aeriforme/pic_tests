import openpmd_api as io
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import micron, centi, pi, c, e, m_e, eV 
import vtk
from vtk.util import numpy_support

MeV = 1e6*eV
lambda_SI = 0.8 * micron # wavelength
omega_SI = 2.0 * pi * c / lambda_SI
n_crit = 1.1 * 1e21 / centi**(3) * (1.e-6/lambda_SI)**(2) #m^-3

path = '/m100_scratch/userexternal/aforment/3DSCANSYNCHR/DLT_a020_l15_d1'



series = io.Series(path+"/diags/Fields_Particles/openpmd_%T.bp", io.Access.read_only)
print("The Series contains {0} iterations:".format(len(series.iterations)))


#for n,i in enumerate(series.iterations,start=0):
for n in (4,):
    i = 3974
    

    print("Iteration number \t {0}".format(i))
    
    j = series.iterations[i]

    print("Current iteration contains {0} meshes:".format(len(j.meshes)))
    for m in j.meshes:
        print("\t {0}".format(m))
    print("")
    print("Current iteration contains {0} particle species:".format(len(j.particles)))
    for ps in j.particles:
        print("\t {0}".format(ps))
        print("With records:")
        for r in j.particles[ps]:
            print("\t {0}".format(r))


    part = j.particles["ele_foam"]
    x = part["position"]["x"]
    y = part["position"]["y"]
    z = part["position"]["z"]
    px = part["momentum"]["x"]
    py = part["momentum"]["y"]
    pz = part["momentum"]["z"]
    w = part["weighting"][io.Mesh_Record_Component.SCALAR]

    x_data = x.load_chunk()
    y_data = y.load_chunk()
    z_data = z.load_chunk()
    px_data = px.load_chunk()
    py_data = py.load_chunk()
    pz_data = pz.load_chunk()   
    w_data = w.load_chunk()

    series.flush()

    
    
    E = m_e*c**2*(np.sqrt(1.+(px_data**2+py_data**2+pz_data**2)/(m_e*c)**2)-1.)/MeV 
    Emax = 1000000 #MeV
    Emin = 10 #MeV
    series.flush()
    index = (E>Emin)&(E<Emax)
    np.savetxt('ELEF_X_Y_Z_PX_PY_PZ_E_t15_%04d.txt' % (i), np.c_[x_data[index]/micron,y_data[index]/micron,z_data[index]/micron,px_data[index]/(m_e*c),py_data[index]/(m_e*c),pz_data[index]/(m_e*c),E[index]],fmt='%.6e')



    
del series 
