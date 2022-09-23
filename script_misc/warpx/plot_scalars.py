import numpy as np
import matplotlib.pyplot as plt

root = '/m100_scratch/userexternal/mgalbiat/3DSCANNICS'
path = root+'/diags/reducedfiles/'


fig,ax = plt.subplots(1,6, figsize=(20,5), dpi=300)

x = np.loadtxt(path + 'FieldEnergy.txt')
ax[0].plot(x[:,1], x[:,2], label='tot')
ax[0].plot(x[:,1], x[:,3], label='E')
ax[0].plot(x[:,1], x[:,4], label='B')
ax[0].set_title('field energy [J]')




x = np.loadtxt(path + 'ParticleEnergy.txt')
ax[1].plot(x[:,1], x[:,3], label='ele foam')
ax[1].plot(x[:,1], x[:,4], label='ion foam')
ax[1].plot(x[:,1], x[:,5], label='pho foam')
ax[1].plot(x[:,1], x[:,6], label='ele subs')
ax[1].plot(x[:,1], x[:,7], label='ion subs')
ax[1].plot(x[:,1], x[:,8], label='pho subs')
ax[1].set_title('partitle energy [J]')
ax[1].set_yscale('log')


x = np.loadtxt(path + 'ParticleNumber.txt')
ax[2].plot(x[:,1], x[:,3], label='ele foam')
ax[2].plot(x[:,1], x[:,4], label='ion foam')
ax[2].plot(x[:,1], x[:,5], label='pho foam')
ax[2].plot(x[:,1], x[:,6], label='ele subs')
ax[2].plot(x[:,1], x[:,7], label='ion subs')
ax[2].plot(x[:,1], x[:,8], label='pho subs')
ax[2].set_title('particle number')
ax[2].set_yscale('log')


x= np.loadtxt(path + 'FieldMaximum.txt')
ax[3].plot(x[:,1], x[:,2], label='Ex')
ax[3].plot(x[:,1], x[:,3], label='Ey')
ax[3].plot(x[:,1], x[:,4], label='Ez')
ax[3].plot(x[:,1], x[:,5], label='|E|')
ax[3].plot(x[:,1], x[:,6], label='Bx')
ax[3].plot(x[:,1], x[:,7], label='By')
ax[3].plot(x[:,1], x[:,8], label='Bz')
ax[3].plot(x[:,1], x[:,9], label='|B|')
ax[3].set_title('field maximum')



x = np.loadtxt(path + 'ParticleExtrema.txt')
ax[4].plot(x[:,1], x[:,18], label='min')
ax[4].plot(x[:,1], x[:,19], label='max')
ax[4].set_title('chi')


ax[5].plot(x[:,1], x[:,14], label='min')
ax[5].plot(x[:,1], x[:,15], label='max')
ax[5].set_title('gamma')




#[0]step() [1]time(s) [2]xmin(m) [3]xmax(m) [4]ymin(m) [5]ymax(m) [6]zmin(m) [7]zmax(m) [8]pxmin(kg*m/s) [9]pxmax(kg*m/s) [10]pymin(kg*m/s) [11]pymax(kg*m/s) [12]pzmin(kg*m/s) [13]pzmax(kg*m/s) [14]gmin() [15]gmax() [16]wmin() [17]wmax() [18]chimin() [19]chimax()


for a in ax.reshape(-1):
    a.legend()
    a.set_xlabel('time [s]')

plt.tight_layout()
plt.savefig('scalars.png')
