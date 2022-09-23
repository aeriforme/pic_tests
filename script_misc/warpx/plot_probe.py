import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
import matplotlib.pylab as pl 
from scipy.constants import c, eV, m_e

MeV = 1e6*eV 
mc2 = m_e*c**2
rc('lines', linewidth=2.0)
rc('font', size=12)

#root = '/m100_scratch/userexternal/mgalbiat/3DSCANNICS'
#path = root+'/diags/reducedfiles/'


fig,ax = plt.subplots(1,1, figsize=(15,15), dpi=300)

path = './'
data = np.loadtxt(path + 'probe.txt')
#cols = pl.cm.rainbow(np.linspace(0,1,steps))

xp = np.unique(data[:,2])
yp = np.unique(data[:,3])
zp = np.unique(data[:,4])

print(yp)
print(zp)

y = data[:,3]
z = data[:,4]
t = data[:,1]
steps = np.shape(data)[0]

S = data[:,10]



#S1 = S[(y==yp[0])&(z==zp[-1])]
#t1 = t[(y==yp[0])&(z==zp[-1])]

for i in yp:
    for j in zp:
        ax.plot(t[(y==i)&(z==j)], S[(y==i)&(z==j)])




#ax[0].set_xscale('log') 
#ax[0].set_yscale('log')

#for t in range(steps):
    

    #ax[0].plot(mc2*bins/MeV, x[t,:],color=cols[t])
#ax[0].set_ylim(1e2,1e15)   
#ax[0].set_yscale('log')
#ax[0].set_title('foam electrons')




#for a in ax.reshape(-1):
#    a.set_xlabel('E [MeV]')
#    a.set_ylabel('dN/dE [arb. units?]')

plt.tight_layout()
plt.savefig('probe.png')
