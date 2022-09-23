import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
import matplotlib.pylab as pl 
from scipy.constants import c, eV, m_e

MeV = 1e6*eV 
mc2 = m_e*c**2
rc('lines', linewidth=2.0)
rc('font', size=12)

root = '/m100_scratch/userexternal/mgalbiat/3DSCANNICS'
path = root+'/diags/reducedfiles/'


fig,ax = plt.subplots(2,3, figsize=(15,5), dpi=300)


nb = 256
bmin = 0
bmax = 400
db = (bmax-bmin)/nb
bins = np.linspace(bmin-db,bmax+db,nb+2)

x = np.loadtxt(path + 'elehistf.txt')
steps = np.shape(x)[0]
cols = pl.cm.rainbow(np.linspace(0,1,steps))
for t in range(steps):
    ax[0][0].plot(mc2*bins/MeV, x[t,:],color=cols[t])
ax[0][0].set_ylim(1e2,1e15)   
ax[0][0].set_yscale('log')
#ax[0][0].set_xscale('log')
ax[0][0].set_title('foam electrons')


x = np.loadtxt(path + 'phohistf.txt')
for t in range(steps):
    ax[1][0].plot(mc2*bins/(MeV), x[t,:],color=cols[t])
#ax[1][0].set_xscale('log')
ax[1][0].set_yscale('log')
ax[1][0].set_title('foam photons')


x = np.loadtxt(path + 'elehists.txt')
for t in range(steps):
    ax[0][1].plot(mc2*bins/(MeV), x[t,:],color=cols[t])
#ax[0][1].set_xscale('log')
ax[0][1].set_yscale('log')
ax[0][1].set_title('subs electrons')



x = np.loadtxt(path + 'phohists.txt')
for t in range(steps):
    ax[1][1].plot(mc2*bins/(MeV), x[t,:],color=cols[t])
#ax[1][1].set_xscale('log')
ax[1][1].set_yscale('log')
ax[1][1].set_title('subs photons')


x = np.loadtxt(path + 'elehistc.txt')
for t in range(steps):
    ax[0][2].plot(mc2*bins/(MeV), x[t,:],color=cols[t])
#ax[0][2].set_xscale('log')
ax[0][2].set_yscale('log')
ax[0][2].set_title('cont electrons')



x = np.loadtxt(path + 'phohistc.txt')
for t in range(steps):
    ax[1][2].plot(mc2*bins/(MeV), x[t,:],color=cols[t])
#ax[1][2].set_xscale('log')
ax[1][2].set_yscale('log')
ax[1][2].set_title('cont photons')




for a in ax.reshape(-1):
    a.set_xlabel('E [MeV]')
    a.set_ylabel('dN/dE [arb. units?]')

plt.tight_layout()
plt.savefig('hist.png')
