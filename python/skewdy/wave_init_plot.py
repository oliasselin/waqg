import os
import numpy as np
import matplotlib.pyplot as plt

scratch_location = '/oasis/scratch/comet/oasselin/temp_project/'
folder = 'niskine/skewdy/'
run = 'initspec/'
location = scratch_location+folder+run

depth = 3000
cm_factor = 100

if not os.path.exists('plots/'+run):
    os.makedirs('plots/'+run)

plt.subplot(1, 2, 1)
u_init    = np.loadtxt(location+'output/init.dat')
plt.plot(u_init[:,1]*cm_factor,u_init[:,0]-depth,label='$\sigma_w = 200 $m')
plt.plot(u_init[:,2]*cm_factor,u_init[:,0]-depth,label='$\sigma_w = 100 $m')
plt.plot(u_init[:,3]*cm_factor,u_init[:,0]-depth,label='$\sigma_w =  50 $m')

plt.title('Initial wave velocity')
plt.xlabel('$u$ (cm/s)')
plt.ylabel('Depth (m)')
plt.grid(color='k', linestyle='-', linewidth=0.1)
plt.ylim(-500,0)


#plt.legend(loc='center right',fontsize='small')
#plt.xlim(0,2e-5)
#plt.ylim(-depth,0)
#plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))


plt.subplot(1, 2, 2)

spec = np.loadtxt(location+'output/init_spec.dat')


plt.loglog(spec[:,0],spec[:,1],label='$\sigma_w = 200 $m')
plt.loglog(spec[:,0],spec[:,2],label='$\sigma_w = 100 $m')
plt.loglog(spec[:,0],spec[:,3],label='$\sigma_w =  50 $m')
plt.title('Vertical energy spectrum')
plt.xlabel('$m$')
plt.ylabel('Energy density')
plt.xlim(1,100)
plt.ylim(1e-12,1e-3)
plt.legend(loc='upper right',fontsize='small')

plt.grid(color='k', linestyle='-', linewidth=0.1)

plt.subplots_adjust(wspace=0.4)
plt.savefig('plots/'+run+'/wave_init.eps',bbox_inches='tight')
#plt.show()
