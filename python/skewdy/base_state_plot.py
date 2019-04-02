import os
import numpy as np
import matplotlib.pyplot as plt

scratch_location = '/oasis/scratch/comet/oasselin/temp_project/'
folder = 'niskine/skewdy/'
run = 'dE60_dt0.1_512/'
location = scratch_location+folder+run

obs_location = 'observations/'
cm_factor = 100  #Set to 100 for converting m to cm

depth = 3000
N0=0.001550529072004
N02=np.square(N0)
U_scale = 0.1*cm_factor

if not os.path.exists('plots/'+run):
    os.makedirs('plots/'+run)

plt.subplot(1, 2, 1)
n2_raw    = np.loadtxt(obs_location+'n2_raw.dat')
n2_smooth = np.loadtxt(obs_location+'n2_smooth.dat')
mid_depth = -np.arange(1.5,2731.5,1)

fit = np.loadtxt(location+'output/rucoeff.dat')
H = depth/(2*np.pi)
full_depth = fit[:,0]*H
full_depth = -full_depth[::-1]
n2_fit = N02*fit[:,2]

plt.plot(n2_raw,mid_depth,'b-',linewidth=.01,label='Raw')
plt.plot(n2_smooth,mid_depth,'k-',linewidth=1.,label='50-day moving average')
plt.plot(n2_fit,full_depth,'r-',linewidth=2.,label='Skewed gaussian fit')
plt.title('Base-state stratification')
plt.xlabel('$N^2$ (s$^{-2}$)')
plt.ylabel('Depth (m)')
plt.grid(color='k', linestyle='-', linewidth=0.1)

plt.legend(loc='center right',fontsize='x-small')
plt.xlim(0,2e-5)
plt.ylim(-depth,0)
plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))


plt.subplot(1, 2, 2)

u_mean = np.loadtxt(location+'output/u_mean.dat')
u_mean = U_scale*u_mean[:,1]
dz_over_2 = fit[0,0]*H
full_depth_s = full_depth-dz_over_2    #Staggered grid (dimensional) 

plt.plot(u_mean,full_depth_s,'r-',linewidth=2.)
plt.title('Base-state velocity')
plt.xlabel('$U$ (cm/s)')
plt.ylabel('Depth (m)')
plt.ylim(-depth,0)

plt.grid(color='k', linestyle='-', linewidth=0.1)

plt.subplots_adjust(wspace=0.4)
plt.savefig('plots/'+run+'/base_state.eps',bbox_inches='tight')
#plt.show()
