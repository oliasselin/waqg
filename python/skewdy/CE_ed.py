#Difference with CE.py: here we have a look at the Eady case. The loss of geostrophic energy is evaluated by subtracting the energy time series of the non-feedback run.

#!/usr/bin/env python                                                                                        
import os
import numpy as np
import matplotlib.pyplot as plt


scratch_location = '/oasis/scratch/comet/oasselin/temp_project/'
folder = 'niskine/skewdy/'
u0='10'
run_fb = 'storm5_uw10/' #'storm5_uw10_ed/'
run_nf = 'storm5/' #'storm5_ed/'
run = 'storm5'

location = scratch_location+folder+run
location_fb = scratch_location+folder+run_fb
location_nf = scratch_location+folder+run_nf

focus_time = 40  #Focus on the first $focus_time days


if not os.path.exists('plots/'+run):
    os.makedirs('plots/'+run)

###################################################
#Vertical profile of the horizontally-averaged WKE#
###################################################


#Load the eddy and coupled energy time series
path_ce = location_fb+'output/ce.dat'
path_ee_fb = location_fb+'output/energy.dat'
path_ee_nf = location_nf+'output/energy.dat'

ce = np.loadtxt(path_ce)
ee_fb = np.loadtxt(path_ee_fb)
ee_nf = np.loadtxt(path_ee_nf)

#Take only the part of the timeseries including focus_time
for index in range(len(ee_fb[:,0])):
    if ee_fb[index,0] > focus_time:
        ee_fb = ee_fb[:index+1,:]
        ee_nf = ee_nf[:index+1,:]
        ce = ce[:index+1,:]
        break

time=ee_fb[:,0]
total_fb = ee_fb[:,1]+ee_fb[:,2]
total_nf = ee_nf[:,1]+ee_nf[:,2]
#ee_perc = total/total[0]

loss = total_nf-total_fb

wpe = ce[:,1]
wce = ce[:,2]
coupled_energy = wpe + wce

norm = 100/total_nf[0]  #Normalize with initial geostrophic (or coupled) energy
#norm=1.


fig = plt.figure(figsize=(8,4))
ax = fig.add_subplot(1, 1, 1)
ax.grid(color='k', linestyle='-', linewidth=0.1)
plt.plot(time,loss*norm,'-k',label='Loss of total geostrophic energy')
plt.plot(time,coupled_energy*norm,'-b',label='WPE+WCE')
plt.plot(time,wpe*norm,'--b',label='WPE')
plt.plot(time,wce*norm,'-.b',label='WCE')
plt.title('Coupled energy conservation, $u_0$ = '+u0+' cm/s')
plt.legend(loc='best',fontsize='small')
plt.xlabel('Time (days)')                                                                                                                                                           
plt.ylabel('Energy (%)')
plt.xlim(0,focus_time)
plt.show()
#plt.savefig('plots/'+run+'/CE.eps',bbox_inches='tight')


