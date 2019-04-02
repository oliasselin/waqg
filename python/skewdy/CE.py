#!/usr/bin/env python                                                                                                                                                               
import os
import numpy as np
import matplotlib.pyplot as plt


scratch_location = '/oasis/scratch/comet/oasselin/temp_project/'
folder = 'niskine/skewdy/'
run = 'storm3_fb_nodiss/'
location = scratch_location+folder+run

focus_time =30  #Focus on the first $focus_time days


if not os.path.exists('plots/'+run):
    os.makedirs('plots/'+run)

###################################################
#Vertical profile of the horizontally-averaged WKE#
###################################################


#Load the eddy energy time series
path_ee = location+'output/energy.dat'
ee = np.loadtxt(path_ee)

#Take only the part of the timeseries including focus_time
for index in range(len(ee[:,0])):
    if ee[index,0] > focus_time:
        ee       =   ee[:index,:]
        break

time=ee[:,0]
total = ee[:,1]+ee[:,2]
ee_perc = total/total[0]

loss = total[0]-total


#Load the coupled energy
path_ce = location+'output/ce.dat'
ce = np.loadtxt(path_ce)
wpe = ce[:,1]
wce = ce[:,2]
coupled_energy = wpe + wce

norm = 100/total[0]  #Normalize with initial geostrophic (or coupled) energy


fig = plt.figure(figsize=(8,4))
ax = fig.add_subplot(1, 1, 1)
ax.grid(color='k', linestyle='-', linewidth=0.1)
plt.plot(time,loss*norm,label='Loss of geostrophic energy')
plt.plot(time,coupled_energy*norm,label='Coupled energy')
plt.plot(time,wpe*norm,label='WPE')
plt.plot(time,wce*norm,label='WCE')
plt.title('Coupled energy conservation, run = '+run)
plt.legend(loc='best')
plt.xlabel('Time (days)')                                                                                                                                                           
plt.ylabel('Energy (%)')
plt.show()
#plt.savefig('plots/'+run+'/EE.eps',bbox_inches='tight')


