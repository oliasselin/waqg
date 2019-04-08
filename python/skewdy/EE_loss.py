#!/usr/bin/env python                                                                                                                                                               
import os
import numpy as np
import matplotlib.pyplot as plt


scratch_location = '/oasis/scratch/comet/oasselin/temp_project/'
folder = 'niskine/skewdy/'

run = 'storm3/'    #Run to compare to
run_list = ['storm3_fb/','storm3_fb_uw2/','storm3_fb_uw4/']
storm_strength = ['10','20','40']


focus_time =10  #Focus on the first $focus_time days


if not os.path.exists('plots/'+run):
    os.makedirs('plots/'+run)


###################################################
#Vertical profile of the horizontally-averaged WKE#
###################################################


#Load the eddy energy time series
path_ee = scratch_location+folder+run+'output/energy.dat'
ee = np.loadtxt(path_ee)

#Take only the part of the timeseries including focus_time
for index in range(len(ee[:,0])):
    if ee[index,0] > focus_time:
        ee   =   ee[:index+1,:]
        break

time=ee[:,0]
control_ee = ee[:,1]+ee[:,2]

loss = np.zeros( (len(time),len(run_list)+1) )
loss[:,0] = time[:]

for run_no,run_fb in enumerate(run_list):

    #Load the eddy energy time series for the feedback run
    path_ee_fb = scratch_location+folder+run_fb+'output/energy.dat'
    ee_fb = np.loadtxt(path_ee_fb)
    ee_fb = ee_fb[:len(time),:]     #Keep only the relevant times
    
    feedback_ee = ee_fb[:,1]+ee_fb[:,2]
        
    loss[:,run_no+1] = 100*(control_ee - feedback_ee)/control_ee[0]  #Percentage of total eddy energy loss


fig = plt.figure(figsize=(8,4))
ax = fig.add_subplot(1, 1, 1)
ax.grid(color='k', linestyle='-', linewidth=0.1)
for run_no,run_fb in enumerate(run_list):
    
    plt.plot(time,loss[:,run_no+1],label='$u_0$ = '+storm_strength[run_no]+' cm/s')

plt.legend(loc='best')
plt.xlabel('Time (days)')                                                                                                                                                           
plt.ylabel('Total geostrophic energy loss (%)')
plt.xlim(0,focus_time)
#plt.show()
plt.savefig('plots/'+run+'/EE_loss.eps',bbox_inches='tight')


