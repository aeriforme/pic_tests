import numpy as np
import matplotlib.pyplot as plt

root = '/m100_scratch/userexternal/aforment/3DSCANSYNCHR/DLT_a020_l15_d1'
path = root+'/diags/reducedfiles/'


fig,ax = plt.subplots(1,1, figsize=(6,5), dpi=300)


x = np.loadtxt(path + 'FieldEnergy.txt')

ax.plot(x[:,1], x[:,2], label='tot')
ax.plot(x[:,1], x[:,3], label='E')
ax.plot(x[:,1], x[:,4], label='B')
ax.set_ylabel('field energy [J]')
ax.set_xlabel('time [s]')

laser_energy=np.max(x[:,2])
ax.set_title('laser energy = %.2f J' % laser_energy)



print('laser energy from PIC = %.2f J' % laser_energy)


from scipy.constants import *


tau = 30 *femto
waist = 3*micron
a0 = 20
I = (a0 / (0.85*0.8))**2*1e18 # W/cm^2
I = I / centi**2  # W/m^2
sigma=pi*waist**2
en = I*tau*sigma

print('laser energy simple stimate = %.2f J' % en)


plt.tight_layout()
plt.savefig('laser_energy.png')
