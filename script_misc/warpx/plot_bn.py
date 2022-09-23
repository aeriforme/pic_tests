import openpmd_api as io
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import Normalize, LogNorm
from scipy.constants import micron, pi, c, centi, eV, e

lambda_SI = 0.8 * micron # wavelength
omega_SI = 2.0 * pi * c / lambda_SI
n_crit = 1.1 * 1e21 / centi**(3) * (1.e-6/lambda_SI)**(2) #m^-3
MeV_SI = 1e6*eV

dx = 1./25.
Ly = 60.

rootm = '/m100_scratch/userexternal/mgalbiat/3DSCANNICS/'
roota = '/m100_scratch/userexternal/aforment/3DSCANSYNCHR/'

#list_dirs=(rootm+'DLT_a020_l2_d1', rootm+'DLT_a020_l5_d1', roota+'DLT_a020_l10_d1', roota+'DLT_a020_l15_d1', roota+'DLT_a020_l20_d1', roota+'DLT_a020_l25_d1') 

list_dirs=(roota+'DLT_a020_l15_d1', roota+'DLT_a020_l20_d1', roota+'DLT_a020_l25_d1')
thick=(15, 20, 25,) 

list_dirs=(rootm+'DLT_a020_l2_d1', rootm+'DLT_a020_l5_d1', roota+'DLT_a020_l10_d1')
thick=(2,5,10)

for k,d in enumerate(list_dirs):
    path=d+'/diags/Fields_Particles'
    print(path)
    series = io.Series(path+"/openpmd_%T.bp", io.Access.read_only)
    #print("The Series contains {0} iterations:".format(len(series.iterations)))

    iter = np.asarray(series.iterations)
    for n,i in enumerate(iter[2:],start=0):
        print("Iteration number \t {0}".format(i))
    
        j = series.iterations[i]
    
        #print("Current iteration contains {0} meshes:".format(len(j.meshes)))
        #for m in j.meshes:
        #    print("\t {0}".format(m))
        #    print("")
        #    print("Current iteration contains {0} particle species:".format(len(j.particles)))
        #    for ps in j.particles:
        #        print("\t {0}".format(ps))

                
        Bz = j.meshes["B"]["z"]
        #By = j.meshes["B"]["y"]
        #Bx = j.meshes["B"]["x"]
        #Ez = j.meshes["E"]["z"]
        #Ey = j.meshes["E"]["y"]
        #Ex = j.meshes["E"]["x"]

        
        rho_ele_foam = j.meshes["rho_ele_foam"][io.Mesh_Record_Component.SCALAR]
        rho_ele_subs = j.meshes["rho_ele_subs"][io.Mesh_Record_Component.SCALAR]

        Bz_data=Bz.load_chunk()
        #By_data=By.load_chunk() #[::2,::2,:]                  
        #Bx_data=Bx.load_chunk() #[::2,::2,:]                
        #Ez_data=Ez.load_chunk() #[::2,::2,:]
        #Ey_data=Ey.load_chunk()
        #Ex_data=Ex.load_chunk()



        rho_ele_foam_data = rho_ele_foam.load_chunk()
        rho_ele_subs_data = rho_ele_subs.load_chunk()
        
        series.flush()
    
        print('shape is: ', np.shape(Bz_data))
        h = int(0.5*np.shape(Bz_data)[0])

        print('plotting...')
        
        my_dpi = 300
        fig, axs = plt.subplots(ncols=2, nrows=1, figsize=(15,3), dpi=300.)
    
        im=axs[0].imshow(Bz_data[h,:,:]/1e3,cmap='bwr',extent=[0,70,-30,30], vmin=-100, vmax=100)
        axs[0].set_title(r'B$_z$')
        plt.colorbar(im, ax=axs[0], orientation='horizontal')
        
        im=axs[1].imshow((rho_ele_subs_data[h,:,:]+rho_ele_foam_data[h,:,:])/(-e*n_crit),cmap='rainbow',extent=[0,70,-30,30], norm=LogNorm(0.1,10))
        axs[1].set_title(r'n$_e$')
        plt.colorbar(im, ax=axs[1], orientation='horizontal')


        #for x in range(np.shape(Bz_data)[2]):
        #    if x*dx < 50: 
        #        endens = Bz_data[h,:,x]**2+Ex_data[h,:,x]**2+Ey_data[h,:,x]**2
        #        index_waist = np.argmin( np.abs(endens - np.max(endens)*0.5) - 1e-6)
               
                #XX = np.append(XX, x*dx)
                #YY = np.append(YY, index_waist*dx-0.5*Ly)
        #        axs[1].scatter(x*dx, index_waist*dx-0.5*Ly, s = 2, zorder = 11, color = 'black' )
                #axs[1].scatter(x*dx, -(index_waist*dx-0.5*Ly), s = 2, zorder = 11, color = 'yellow' ) #SIMMETRIZZATO

                
        #axs[1].plot(XX,YY, lw=2, zorder=11, color='black') 
                

        axs[0].set_xlim(20,70)
        axs[1].set_xlim(20,70)
        axs[0].set_ylim(-5,5)
        axs[1].set_ylim(-5,5)
#        fig.suptitle('t = {0}'.format(j.time))
        
        img_file_name ='./EB_%02d_%06d.png' % (thick[k], i) 
        plt.tight_layout() #pad=2., h_pad=2., w_pad=2.)
        fig.suptitle('t = {0}'.format(j.time))

        plt.savefig(img_file_name, dpi=300)
        
        plt.close('all')
    
    del series 
