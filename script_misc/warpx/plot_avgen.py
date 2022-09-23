import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
from scipy.constants import femto, eV
import matplotlib.pylab as pl

MeV = 1e6*eV 

rc('lines', linewidth=4.0)
rc('font', size=12)

##[0]step() [1]time(s) [2]total(J) [3]ele_foam(J) [4]ion_foam(J) [5]pho_foam(J) [6]ele_subs(J) [7]ion_subs(J) [8]pho_subs(J) [9]ele_cont(J) [10]ion_cont(J) [11]pho_cont(J) [12]total_mean(J) [13]ele_foam_mean(J) [14]ion_foam_mean(J) [15]pho_foam_mean(J) [16]ele_subs_mean(J) [17]ion_subs_mean(J) [18]pho_subs_mean(J) [19]ele_cont_mean(J) [20]ion_cont_mean(J) [21]pho_cont_mean(J)

rootm = '/m100_scratch/userexternal/mgalbiat/3DSCANNICS/'
roota = '/m100_scratch/userexternal/aforment/3DSCANSYNCHR/'

dirs = (rootm+'DLT_a020_l0_d1', rootm+'DLT_a020_l2_d1', rootm+'DLT_a020_l5_d1', roota+'DLT_a020_l10_d1', roota+'DLT_a020_l15_d1', roota+'DLT_a020_l20_d1', roota+'DLT_a020_l25_d1') 
labs = ('0', '2', '5', '10', '15', '20', '25') 
fig,ax = plt.subplots(3,3, figsize=(15,15), dpi=300)
cols = pl.cm.rainbow(np.linspace(0,1,len(dirs)))


for i,d in enumerate(dirs):

    path =  d + '/diags/reducedfiles/'

    x = np.loadtxt(path + 'ParticleEnergy.txt')


    if i == 0 :

        # #[0]step() [1]time(s) [2]total(J) [3]ele_subs(J) [4]ion_subs(J) [5]pho_subs(J) [6]ele_cont(J) [7]ion_cont(J) [8]pho_cont(J) [9]total_mean(J) [10]ele_subs_mean(J) [11]ion_subs_mean(J) [12]pho_subs_mean(J) [13]ele_cont_mean(J) [14]ion_cont_mean(J) [15]pho_cont_mean(J)

        ax[0][1].plot(x[:,1]/femto, x[:,10]/MeV, label=labs[i], color=cols[i])
        ax[0][1].set_title('ele subs')

        ax[0][2].plot(x[:,1]/femto, x[:,13]/MeV, label=labs[i], color=cols[i])
        ax[0][2].set_title('ele cont')

        ax[1][1].plot(x[:,1]/femto, x[:,11]/MeV, label=labs[i], color=cols[i])
        ax[1][1].set_title('ion subs')

        ax[1][2].plot(x[:,1]/femto, x[:,14]/MeV, label=labs[i], color=cols[i])
        ax[1][2].set_title('ion cont')

        ax[2][1].plot(x[:,1]/femto, x[:,12]/MeV, label=labs[i], color=cols[i])
        ax[2][1].set_title('pho subs')

        ax[2][2].plot(x[:,1]/femto, x[:,15]/MeV, label=labs[i], color=cols[i])
        ax[2][2].set_title('pho cont')
        
        continue 

    ax[0][0].plot(x[:,1]/femto, x[:,13]/MeV, label=labs[i], color=cols[i])
    ax[0][0].set_title('ele foam')

    ax[0][1].plot(x[:,1]/femto, x[:,16]/MeV, label=labs[i], color=cols[i])
    ax[0][1].set_title('ele subs')

    ax[0][2].plot(x[:,1]/femto, x[:,19]/MeV, label=labs[i], color=cols[i])
    ax[0][2].set_title('ele cont')

    ax[1][0].plot(x[:,1]/femto, x[:,14]/MeV, label=labs[i], color=cols[i])
    ax[1][0].set_title('ion foam')
    
    ax[1][1].plot(x[:,1]/femto, x[:,17]/MeV, label=labs[i], color=cols[i])
    ax[1][1].set_title('ion subs')

    ax[1][2].plot(x[:,1]/femto, x[:,20]/MeV, label=labs[i], color=cols[i])
    ax[1][2].set_title('ion cont')

    ax[2][0].plot(x[:,1]/femto, x[:,15]/MeV, label=labs[i], color=cols[i])
    ax[2][0].set_title('pho foam')

    ax[2][1].plot(x[:,1]/femto, x[:,18]/MeV, label=labs[i], color=cols[i])
    ax[2][1].set_title('pho subs')

    ax[2][2].plot(x[:,1]/femto, x[:,21]/MeV, label=labs[i], color=cols[i])
    ax[2][2].set_title('pho cont')

    
    for a in ax.reshape(-1):
        a.set_xlabel('time [fs]')
        a.set_ylabel('species average energy [MeV]')

    ax[1][2].legend(title=r't [$\mu$m]')

    plt.tight_layout()
    plt.savefig('avgen.png')
