import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc

rc('lines', linewidth=4.0)
rc('font', size=12)


rootm = '/m100_scratch/userexternal/mgalbiat/3DSCANNICS/'
roota = '/m100_scratch/userexternal/aforment/3DSCANSYNCHR/'

dirs = (rootm+'DLT_a020_l0_d1', rootm+'DLT_a020_l2_d1', rootm+'DLT_a020_l5_d1', roota+'DLT_a020_l10_d1', roota+'DLT_a020_l15_d1', roota+'DLT_a020_l20_d1', roota+'DLT_a020_l25_d1') 
labs = ('0', '2', '5', '10', '15', '20', '25') 

fig,ax = plt.subplots(1,3, figsize=(15,5), dpi=300)



cols = ('#003f5c', '#374c80', '#7a5195','#bc5090','#ef5675','#ff764a','#ffa600')

for i,d in enumerate(dirs):

    path =  d + '/diags/reducedfiles/'

    x = np.loadtxt(path + 'FieldEnergy.txt')
    laser_energy =  np.max(x[:,2])

    x = np.loadtxt(path + 'ParticleEnergy.txt')
    ax[0].plot(x[:,1], (x[:,3]+x[:,6]+x[:,9])/laser_energy, label=labs[i], color=cols[i])
    ax[0].set_title('electrons')


    ax[1].plot(x[:,1], (x[:,4]+x[:,7]+x[:,10])/laser_energy, label=labs[i], color=cols[i])
    ax[1].set_title('ions')

    ax[2].plot(x[:,1], (x[:,5]+x[:,8]+x[:,11])/laser_energy, label=labs[i], color=cols[i])
    ax[2].set_title('photons')

    for a in ax.reshape(-1):
        a.legend(title=r't [$\mu$m]')
        a.set_xlabel('time [s]')
        a.set_ylabel('particle energy / laser energy')
        
    plt.tight_layout()
    plt.savefig('conveff.png')
