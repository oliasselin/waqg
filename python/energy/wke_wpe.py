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
path_energy=location+'/output/we.dat'
if os.path.isfile(path_energy) and os.stat(path_energy).st_size != 0:
    we_nofb=np.loadtxt(path_energy)                                      #0: time (days)      1: wke total   2: wpe total   3: wke in kh=0 mode                          

run = 'storm_Uw0.1_fb'
location = scratch_location+folder+run
path_energy=location+'/output/we.dat'
if os.path.isfile(path_energy) and os.stat(path_energy).st_size != 0:
    we_fb10=np.loadtxt(path_energy)                                      #0: time (days)      1: wke total   2: wpe total   3: wke in kh=0 mode                          

run = 'storm_Uw0.2_fb'
location = scratch_location+folder+run
path_energy=location+'/output/we.dat'
if os.path.isfile(path_energy) and os.stat(path_energy).st_size != 0:
    we_fb20=np.loadtxt(path_energy)                                      #0: time (days)      1: wke total   2: wpe total   3: wke in kh=0 mode                          


#Normalize energy with the initial value of wave kinetic energy: 
kinit = we_nofb[0,1]
we_nofb[:,2]=we_nofb[:,2]/kinit
we_nofb[:,1]=we_nofb[:,1]/kinit

kinit = we_fb10[0,1]
we_fb10[:,2]=we_fb10[:,2]/kinit
we_fb10[:,1]=we_fb10[:,1]/kinit

kinit = we_fb20[0,1]
we_fb20[:,2]=we_fb20[:,2]/kinit
we_fb20[:,1]=we_fb20[:,1]/kinit




fig, ax = plt.subplots(1)
plt.plot(we_nofb[:,0],we_nofb[:,1],color='b',linestyle='-' ,label='WKE, no feedback')
plt.plot(we_nofb[:,0],we_nofb[:,2],color='b',linestyle='--',label='WPE, no feedback')
plt.plot(we_fb10[:,0],we_fb10[:,1],color='g',linestyle='-', label='WKE, $U_w = 10$ cm/s')
plt.plot(we_fb10[:,0],we_fb10[:,2],color='g',linestyle='--',label='WPE, $U_w = 10$ cm/s')
plt.plot(we_fb20[:,0],we_fb20[:,1],color='r',linestyle='-', label='WKE, $U_w = 20$ cm/s')
plt.plot(we_fb20[:,0],we_fb20[:,2],color='r',linestyle='--',label='WPE, $U_w = 20$ cm/s')
ax.set_xlabel('Time (days)')
ax.set_ylabel('Normalized Energy')
ax.set_xlim([0,30])
plt.legend(prop={'size': 10})
plt.title('Wave kinetic and potential energies')
fig.savefig('plots/wke_wpe.png')
plt.close(fig)
#plt.show()     
