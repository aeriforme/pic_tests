import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
import matplotlib.pylab as pl 
from scipy.constants import c, eV, m_e

MeV = 1e6*eV 
mc2 = m_e*c**2
rc('lines', linewidth=2.0)
rc('font', size=12)

rootm = '/m100_scratch/userexternal/mgalbiat/3DSCANNICS/'
roota = '/m100_scratch/userexternal/aforment/3DSCANSYNCHR/'

dirs = (rootm+'DLT_a020_l0_d1', rootm+'DLT_a020_l2_d1', rootm+'DLT_a020_l5_d1', roota+'DLT_a020_l10_d1', roota+'DLT_a020_l15_d1', roota+'DLT_a020_l20_d1', roota+'DLT_a020_l25_d1') 
labs = ('0', '2', '5', '10', '15', '20', '25') 


fig,ax = plt.subplots(2,3, figsize=(15,5), dpi=300)


nb = 512
bmin = 0
bmax = 400
db = (bmax-bmin)/nb
bins = np.linspace(bmin-db,bmax+db,nb+2)
steps = 153
t = steps-1
cols = pl.cm.rainbow(np.linspace(0,1,len(dirs)))

for i,d in enumerate(dirs):

    path =  d + '/diags/reducedfiles/'

    if i > 0: 
        x = np.loadtxt(path + 'elehistf.txt')
        ax[0][0].plot(mc2*bins/MeV, x[t,:],color=cols[i])
        ax[0][0].set_ylim(1e2,1e15)   
        ax[0][0].set_yscale('log')
        #ax[0][0].set_xscale('log')
        ax[0][0].set_title('foam electrons')


        x = np.loadtxt(path + 'phohistf.txt')
        ax[1][0].plot(mc2*bins/(MeV), x[t,:],color=cols[i])
        #ax[1][0].set_xscale('log')
        ax[1][0].set_yscale('log')
        ax[1][0].set_title('foam photons')


    x = np.loadtxt(path + 'elehists.txt')
    ax[0][1].plot(mc2*bins/(MeV), x[t,:],color=cols[i])
    #ax[0][1].set_xscale('log')
    ax[0][1].set_yscale('log')
    ax[0][1].set_title('subs electrons')

    x = np.loadtxt(path + 'phohists.txt')
    ax[1][1].plot(mc2*bins/(MeV), x[t,:],color=cols[i])
    #ax[1][1].set_xscale('log')
    ax[1][1].set_yscale('log')
    ax[1][1].set_title('subs photons')

    x = np.loadtxt(path + 'elehistc.txt')
    ax[0][2].plot(mc2*bins/(MeV), x[t,:],color=cols[i])
    #ax[0][2].set_xscale('log')
    ax[0][2].set_yscale('log')
    ax[0][2].set_title('cont electrons')


    x = np.loadtxt(path + 'phohistc.txt')
    ax[1][2].plot(mc2*bins/(MeV), x[t,:],color=cols[i], label=labs[i])
    #ax[1][2].set_xscale('log')
    ax[1][2].set_yscale('log')
    ax[1][2].set_title('cont photons')

for a in ax.reshape(-1):
    a.set_xlabel('E [MeV]')
    a.set_ylabel('dN/dE [arb. units?]')

plt.tight_layout()
plt.savefig('hist.png')
