#!/usr/bin/env python                                                                                                                                               
import os
import numpy as np
import matplotlib.pyplot as plt



#Location of things on the machine                                                                                                                                  
scratch_location = '/oasis/scratch/comet/oasselin/temp_project/'
folder = 'niskine/expeady/'
plot_dir = 'plots/'
if not os.path.exists(plot_dir):
    os.makedirs(plot_dir)




run = 'storm_Uw0.2_2/'
location = scratch_location+folder+run
path_energy=location+'/output/energy.dat'
if os.path.isfile(path_energy) and os.stat(path_energy).st_size != 0:
    eke_nofb=np.loadtxt(path_energy)                                      #0: time (days)      1: wke total   2: wpe total   3: wke in kh=0 mode                          

run = 'storm_Uw0.2_fb2'
location = scratch_location+folder+run
path_energy=location+'/output/energy.dat'
if os.path.isfile(path_energy) and os.stat(path_energy).st_size != 0:
    eke_fb20=np.loadtxt(path_energy)                                      #0: time (days)      1: wke total   2: wpe total   3: wke in kh=0 mode                          

path_energy=location+'/output/we.dat'
if os.path.isfile(path_energy) and os.stat(path_energy).st_size != 0:
    we_fb20=np.loadtxt(path_energy)                                      #0: time (days)      1: wke total   2: wpe total   3: wke in kh=0 mode                          



fig, ax = plt.subplots(1)
plt.plot(eke_nofb[:,0],eke_nofb[:,1],color='k',linestyle='-', label='EKE, no feedback')
plt.plot(eke_fb20[:,0],eke_fb20[:,1],color='c',linestyle='-',label='EKE, $U_w = 20$ cm/s')
plt.plot(we_fb20[:,0],we_fb20[:,1],color='r',linestyle='-', label='WKE, $U_w = 20$ cm/s')
plt.plot(we_fb20[:,0],we_fb20[:,2],color='r',linestyle='--',label='WPE, $U_w = 20$ cm/s')
ax.set_xlabel('Time (days)')
ax.set_ylabel('Energy (m/s)$^2$')
ax.set_xlim([0,30])
plt.legend(prop={'size': 10},loc=7)
plt.title('Eddy and wave energies')
fig.savefig('plots/eke_we.png')
plt.close(fig)
#plt.show()     
