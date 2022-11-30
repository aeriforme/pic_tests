import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import micron,c,pi,centi,femto,e,epsilon_0,m_e
import openpmd_api as io
import happi 
import sdf 
from matplotlib import use

#use('AGG')

#__________________________________________________________
# create plot directory
plot_dir = './plots'
if os.path.exists(plot_dir) is False:
    os.mkdir(plot_dir)


#__________________________________________________________
lambda_SI = 0.8*micron # wavelength
omega_SI = 2.0*pi*c / lambda_SI
fs = 1.e-15 * omega_SI;
mc2 = 0.510998950e6 
my_dpi = 300

#__________________________________________________________
# simulation parameters
Lx = 70*micron
resx = 20. # points per micron 
nx = Lx * resx / micron 
dx = Lx/nx
micron_s = 1.e-6 * omega_SI/c
vol = 70*30*micron_s**(2)
n_crit = (m_e*epsilon_0*(2*pi*c)**2)/((e*lambda_SI)**2)
nppc_ic = 16.
delta_E = 40/0.510998950/200
const = delta_E*vol*n_crit*(c/omega_SI)**2

cfl = 0.98
dt = cfl * dx / (c * np.sqrt(2.))

Tsim = 1.5*Lx / c 
nsteps = int(Tsim/dt) 


every_fs = np.floor(femto/dt)
out_freq = int(10*every_fs)

steps = out_freq * np.arange(nsteps)
times = steps * dt/femto

steps = np.arange(0,nsteps,out_freq)

#warpx
warpx_dir='./warpx/diags/reducedfiles'
data_w = np.loadtxt(warpx_dir+'/FieldEnergy.txt')  
times_w = data_w[:,1]/femto
Uelm_w = data_w[:,2]

Ukin_w = np.loadtxt(warpx_dir+'/ParticleEnergy.txt')
Ukin_ef_w = Ukin_w[:,3]
Ukin_es_w = Ukin_w[:,4]
Ukin_ic_w = Ukin_w[:,8]

datah_w = np.loadtxt(warpx_dir+'/ParticleHist.txt')[1:,1:]
ene_w = np.linspace(0,40,200)
t_w = datah_w[:,1]

#smilei
smilei_dir = './smilei'
s = happi.Open(smilei_dir)
den_s = s.ParticleBinning(2,units=['Mev','um','fs'])
t_s = den_s.getAvailableTimesteps()

#epoch
epoch_dir = './epoch/Data/'
#sdf.list_variables(data)
filenames=np.genfromtxt(epoch_dir+'histogram.visit',dtype='str')
t_e=len(filenames)


for T,T_s in zip(range(np.size(t_w)),t_s):
    #warpx
    den_w = datah_w[T,1:]
    #smilei
    den_s = s.ParticleBinning(diagNumber="#2",timesteps=T_s,units=['um','fs'])
    den_s = den_s.getData()
    #epoch
    fname=epoch_dir+'hist%04d.sdf' % T
    raw=sdf.read(fname)
    den_e = raw.__dict__["dist_fn_enehist_ion_cont"].data

    my_dpi = 300 
    fig, ax = plt.subplots(ncols=1, nrows=1, figsize=(4000./my_dpi, 1000./my_dpi), dpi=my_dpi, sharex=True)
    ax.plot(ene_w,den_w,label='Warpx')
    ax.plot(ene_w,np.transpose(np.array(den_s))*const,label='Smilei')
    ax.plot(ene_w,den_e,label='Epoch')
    ax.set_title('Energy density')
    ax.set_xlabel('Energy [Mev]')
    ax.set_ylabel('Particle density')
    image_file_name =plot_dir+'/energy_den_%04d.png' %T 
    plt.tight_layout()
    plt.savefig(image_file_name,dpi=my_dpi)
    plt.close()

##smilei 
data_s = s.Scalar(scalar='Uelm', units=['J/m', 'fs'])
times_s = (data_s.getTimes())
Uelm_s = data_s.getData()

Ukin_ef_s = s.Scalar(scalar='Ukin_ele_f', units=['J/m', 'fs']).getData()
Ukin_es_s = s.Scalar(scalar='Ukin_ele_s', units=['J/m', 'fs']).getData()
Ukin_ic_s = s.Scalar(scalar='Ukin_ion_c', units=['J/m', 'fs']).getData()

# epoch
Uelm_e=[]
times_e=[]
Ukin_ef_e = []
Ukin_es_e = []
Ukin_ic_e = []


#sh.list_variables(data)
filenames=np.genfromtxt(epoch_dir+'scalars.visit',dtype='str')
ndumps=len(filenames) 

for t in range(ndumps):
    fname=epoch_dir+'scalars%04d.sdf' % t
    data=sdf.read(fname)
    times_e = np.append(times_e,data.Header['time']/femto)
    Uelm_e=np.append(Uelm_e,data.Total_Field_Energy_in_Simulation__J_.data)     
    Ukin_ef_e=np.append(Ukin_ef_e,data.Total_Particle_Energy_ele_foam__J_.data)
    Ukin_es_e=np.append(Ukin_es_e,data.Total_Particle_Energy_ele_subs__J_.data)     
    Ukin_ic_e=np.append(Ukin_ic_e,data.Total_Particle_Energy_ion_cont__J_.data)     
#__________________________________________________________
# plot images
my_dpi = 300 
fig, ax = plt.subplots(ncols=4, nrows=1, figsize=(4000./my_dpi, 1000./my_dpi), dpi=my_dpi, sharex=True)

ax[0].plot(times_w, Uelm_w, label='WarpX')
ax[0].plot(times_s, Uelm_s, label='Smilei')
#ax[0].plot(times_e, Uelm_e, label='EPOCH')
ax[0].set_title('field energy')

ax[1].plot(times_w, Ukin_ef_w, label='WarpX')
ax[1].plot(times_s, Ukin_ef_s, label='Smilei')
#ax[1].plot(times_e, Ukin_ef_e, label='EPOCH')
ax[1].set_title('foam electron energy')

ax[2].plot(times_w, Ukin_es_w, label='WarpX')
ax[2].plot(times_s, Ukin_es_s, label='Smilei')
#ax[2].plot(times_e, Ukin_es_e, label='EPOCH')
ax[2].set_title('subs electron energy')

ax[3].plot(times_w, Ukin_ic_w, label='WarpX')
ax[3].plot(times_s, Ukin_ic_s, label='Smilei')
#ax[3].plot(times_e, Ukin_ic_e, label='EPOCH')
ax[3].set_title('proton energy')

for a in ax.reshape(-1):
    a.set_xlabel('time [fs]')
    a.set_ylabel('energy [J/m]')
    a.legend()

image_file_name =plot_dir+'/scalars.png' 
plt.tight_layout()
plt.savefig(image_file_name,dpi=my_dpi)
plt.close()


