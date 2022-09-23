import openpmd_api as io
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import Normalize
import vtk 
from vtk.util import numpy_support
from scipy.constants import * 


root = '/m100_scratch/userexternal/aforment/3DSCANSYNCHR/DLT_a020_l15_d1'
path = root+'/diags/Fields_Particles'


lambda_SI = 0.8 * micron # wavelength
omega_SI = 2.0 * pi * c / lambda_SI
n_crit = 1.1 * 1e21 / centi**(3) * (1.e-6/lambda_SI)**(2) #m^-3


dx=1./25.
dy=dx
dz=dy

ox=0.
oy=-30.
oz=-30.

series = io.Series(path+"/openpmd_%T.bp", io.Access.read_only)
print("The Series contains {0} iterations:".format(len(series.iterations)))


def write(data,outfile,dx,dy,dz,ox,oy,oz):
    data = np.asarray(data)
    imdata = vtk.vtkImageData()
    depthArray = numpy_support.numpy_to_vtk(np.ravel(data, order='F'),  deep=True, array_type=vtk.VTK_DOUBLE)
    imdata.SetDimensions(data.shape)
    imdata.SetSpacing([dx,dy,dz])
    imdata.SetOrigin([ox,oy,oz])
    imdata.GetPointData().SetScalars(depthArray)
    writer = vtk.vtkMetaImageWriter()
    writer.SetFileName(outfile)
    writer.SetInputData(imdata)
    writer.Write()


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
        #print("With records:")
        #for r in j.particles[ps]:
            #print("\t {0}".format(r))


    Bz = j.meshes["B"]["z"]
    Bz_data=Bz.load_chunk()
    series.flush()

    outfile="bz_t15_%05d.mhd" % i
    data = Bz_data[::3,::3,::3]/1e3
    write(np.transpose(data),outfile,dx*3,dy*3,dz*3,ox,oy,oz)   

    del Bz, Bz_data, data 
    
    rho_ele_foam = j.meshes["rho_ele_foam"][io.Mesh_Record_Component.SCALAR]   
    rho_ele_subs = j.meshes["rho_ele_subs"][io.Mesh_Record_Component.SCALAR]
    rho_ele_foam_data = rho_ele_foam.load_chunk()
    rho_ele_subs_data = rho_ele_subs.load_chunk()
    series.flush()

    outfile="ne_t15_%05d.mhd" % i
    data = -(rho_ele_foam_data[::3,::3,::3]+rho_ele_subs_data[::3,::3,::3])/(e*n_crit)

    write(np.transpose(data),outfile,dx*3,dy*3,dz*3,ox,oy,oz)
    
    
del series 
