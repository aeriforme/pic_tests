import openpmd_api as io
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import Normalize

path = './diags/Fields_Particles'


series = io.Series(path+"/openpmd_%T.bp", io.Access.read_only)
print("The Series contains {0} iterations:".format(len(series.iterations)))


for n,i in enumerate(series.iterations[2:],start=0):
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


    #Ex = j.meshes["E"]["x"]
    #Ey = j.meshes["E"]["y"]
    #Ez = j.meshes["E"]["z"]
    #Bx = j.meshes["B"]["x"]
    #By = j.meshes["B"]["y"]
    #Bz = j.meshes["B"]["z"]
    rho_ele_foam = j.meshes["rho_ele_foam"][io.Mesh_Record_Component.SCALAR]
    rho_ele_subs = j.meshes["rho_ele_subs"][io.Mesh_Record_Component.SCALAR]

    rho_ion_foam = j.meshes["rho_ion_foam"][io.Mesh_Record_Component.SCALAR]
    rho_ion_subs = j.meshes["rho_ion_subs"][io.Mesh_Record_Component.SCALAR]


    #Ex_data=Ex.load_chunk()
    #Ey_data=Ey.load_chunk()
    #Ez_data=Ez.load_chunk()
    #Bx_data=Bx.load_chunk()
    #By_data=By.load_chunk()
    #Bz_data=Bz.load_chunk()

    rho_ele_foam_data = rho_ele_foam.load_chunk()
    rho_ele_subs_data = rho_ele_subs.load_chunk()

    rho_ion_foam_data = rho_ion_foam.load_chunk()
    rho_ion_subs_data = rho_ion_subs.load_chunk()

    
    series.flush()
    

    h = int(0.5*np.shape(rho_ele_foam_data)[0])
    my_dpi = 300
    fig, axs = plt.subplots(ncols=6, nrows=1)
    fig.set_size_inches(6000./my_dpi, 1000./my_dpi)
    fig.set_dpi(my_dpi)

    
    axs[2].imshow(rho_ele_foam_data[h,:,:],cmap='viridis')
    axs[3].imshow(rho_ele_subs_data[h,:,:],cmap='viridis')
    axs[0].imshow(rho_ion_foam_data[h,:,:],cmap='viridis')
    axs[1].imshow(rho_ion_subs_data[h,:,:],cmap='viridis')

    fig.suptitle('t = {0}'.format(j.time))
    
    img_file_name ='./rho_%06d.png' % i 
    plt.tight_layout(pad=2., h_pad=2., w_pad=2.)
    plt.savefig(img_file_name, dpi=300)

    plt.close('all')
    
del series 
