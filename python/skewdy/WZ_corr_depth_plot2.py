import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid
import numpy as np
import os
from finds import find_resolution
from finds import find_scales
from finds import find_timestep

focus_time = 40

run_list=['storm5_uw10/','storm7_uw10/','storm7_uw10/','storm7_uw10/',]
sliceloc=['htop','htop','hmid','hbot',]
depth_list=['Surface','50 m','100 m','200 m']

colors_list=['k','b','g','r']
sn_percentile = 20
    
fig = plt.figure(figsize=(8,4))
ax = fig.add_subplot(1, 1, 1)
ax.grid(color='k', linestyle='-', linewidth=0.1)

for run_no,run_iter in enumerate(run_list):



    out = np.loadtxt('data/'+run_iter+'/WKE_zeta_corr_'+sliceloc[run_no]+'.dat')

    no_pixels=512*512

    plt.plot(out[1:,0],out[1:,4]/no_pixels,'-'+colors_list[run_no],label=depth_list[run_no])
#    plt.plot(out[1:,0],out[1:,3],'-'+color_list[run_iter],label='Full average, '+depth[run_iter])
#    plt.plot(out[1:,0],out[1:,4],'--'+color_list[run_iter],label='only strongly negative $\zeta$, '+depth[run_iter])
    
#plt.plot(out[:,0],out[:,0]*0,'-k',linewidth=0.1,label=None)
plt.legend(loc='best',fontsize='small')
plt.xlabel('Time (days)')
#plt.ylabel('Correlation')
plt.ylabel(r'$\overline{WKE}^{ SN}$')
#plt.ylabel('Average WKE in SN regions (m$^2$/s$^2$)')
plt.xlim(0,focus_time)
#plt.ylim(-1,2)
#plt.show()                                                                                                                                                                           
plt.savefig('plots/'+run_list[0]+'/WKE_depth.eps',bbox_inches='tight')
