#!/usr/bin/env python                                                                                                                                                               
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker
from finds import find_resolution
from finds import find_scales
from finds import find_timestep

show=0
scratch_location = '/oasis/scratch/comet/oasselin/temp_project/'
folder = 'niskine/skewdy/'
run = 'storm5/'
run_spec = 'khspec/' #khspec#
location_spec = scratch_location+folder+run_spec

depth_list  = ['1500','1000','500','400','200','100','50','25','10','0']     #Depths of spectra in meters
which_depth = ['9','7','6','5','4','3','2','1']                                            #Which depths to plot

lines = ["-","--","-.",":","-","--","-.",":"]

colormap = plt.cm.get_cmap('RdBu_r')
#color_sn = colormap(0)
#color_sp = colormap(255)
colors = [colormap(0),colormap(16),colormap(32),colormap(48),colormap(64),colormap(80),colormap(96),colormap(112)]



ymax = 10.
ymin = 1e-5

#Read parameters from the source#
n1,n2,n3 = find_resolution(location_spec)


if not os.path.exists('plots/'+run):
    os.makedirs('plots/'+run)

fig = plt.figure(figsize=(6,6))
ax1 = fig.add_subplot(1, 1, 1)
ax1.grid(color='k', linestyle='-', linewidth=0.1)


for no,depth_no in enumerate(which_depth):
    
    path_spec = location_spec+'output/h'+depth_no+'.dat'

    spec = np.loadtxt(path_spec)

    spec=spec[0:int(int(n1)/3)+1,:]     #Take the first time step for the spetrum

#    ax1.loglog(spec[:,0],spec[:,1],lines[no]+'b',label=depth_list[int(depth_no)]+' m')
    ax1.loglog(spec[:,0],spec[:,1],color=colors[no],label=depth_list[int(depth_no)]+' m')


ax1.loglog(spec[:,0],100*np.power(spec[:,0],-2.),'-k',label='$k_h^{-2}$')

ax1.set_xlabel('$k_h$')
ax1.set_ylabel('Kinetic Energy Density (m$^3$/s$^2$)')
ax1.set_title('Initial EKE spectrum')
ax1.legend(loc='best',fontsize='small')
ax1.set_xlim((1,int(int(n3)/3)))                                                                                                                                             
ax1.set_ylim((ymin,ymax))
if (show==1):
    plt.show()
else:
    plt.savefig('plots/'+run+'/init_khspec2.eps',bbox_inches='tight')



