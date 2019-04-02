#!/usr/bin/env python                                                                                                                                                               
import os
import numpy as np
import matplotlib.pyplot as plt


scratch_location = '/oasis/scratch/comet/oasselin/temp_project/'
folder = 'niskine/skewdy/'
run = 'storm0/'
run_fb = 'storm0_fb/'
location = scratch_location+folder+run
location_fb = scratch_location+folder+run_fb

focus_time =100  #Focus on the first $focus_time days


if not os.path.exists('plots/'+run):
    os.makedirs('plots/'+run)

###################################################
#Vertical profile of the horizontally-averaged WKE#
###################################################


#Load the eddy energy time series
path_ee = location+'output/energy.dat'
ee = np.loadtxt(path_ee)

path_ee = location_fb+'output/energy.dat'
ee_fb = np.loadtxt(path_ee)

#Take only the part of the timeseries including focus_time
for index in range(len(ee[:,0])):
    if ee[index,0] > focus_time:
        ee   =   ee[:index,:]
        ee_fb=ee_fb[:index,:]
        break

time=ee[:,0]
total = ee[:,1]+ee[:,2]
ee_perc = total/total[0]

total_fb = ee_fb[:,1]+ee_fb[:,2]
ee_perc_fb = total_fb/total_fb[0]


fig = plt.figure(figsize=(8,4))
ax = fig.add_subplot(1, 1, 1)
ax.grid(color='k', linestyle='-', linewidth=0.1)
plt.plot(time,ee_perc   ,label='No feedback')
plt.plot(time,ee_perc_fb,label='w/ feedback')
plt.legend()
plt.xlabel('Time (days)')                                                                                                                                                           
plt.ylabel('Geostrophic energy (%)')
#plt.show()
plt.savefig('plots/'+run+'/EE.eps',bbox_inches='tight')


