#!/usr/bin/env python                                                                                                                                                               
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker
from finds import find_resolution
from finds import find_scales
from finds import find_timestep

just_power=1    #0: Plot power and total energy geostrophic loss. 1: Just power!
show=0

scratch_location = '/oasis/scratch/comet/oasselin/temp_project/'
folder = 'niskine/skewdy/'

run = 'storm5/'    #Run to compare to                                                                                                                                                   
location = scratch_location+folder+run
 
run_list = ['storm5_uw10/','storm5_uw20/','storm5_uw40/']
storm_strength = ['10','20','40']
linestyle_list = ['-','-.',':']

focus_time =10  #Focus on the first $focus_time days  



#Read parameters from the source#
n1,n2,n3 = find_resolution(location)
Dx,Dz,L_scale,H_scale,U_scale,h_thermo = find_scales(location)
delt,freq_etot,freq_ez,freq_we,freq_wz = find_timestep(location)

T_scale = L_scale/U_scale
s_to_day=1./(3600*24)
ts_ez_days = delt*freq_ez*T_scale*s_to_day
max_time_ez = np.ceil(focus_time/ts_ez_days)+1#number of time steps to reach focus_time days...
interest_depth = n3-1   #We only care about the surface value


if not os.path.exists('plots/'+run):
    os.makedirs('plots/'+run)


#####################
# Total energy loss #
#####################

#Load the eddy energy time series                                                                                                                                                        
path_ee = scratch_location+folder+run+'output/energy.dat'
ee = np.loadtxt(path_ee)

#Take only the part of the timeseries including focus_time                                                                                                                               
for index in range(len(ee[:,0])):
    if ee[index,0] > focus_time:
        ee   =   ee[:index+1,:]
        break

time=ee[:,0]
dt = time[1]
control_ee = ee[:,1]+ee[:,2]

loss = np.zeros( (len(time),len(run_list)+1) )
loss[:,0] = time[:]
power = np.zeros( (len(time)-1,len(run_list)+1) )
power[:,0] = time[:-1]+dt/2.

for run_no,run_fb in enumerate(run_list):

    #Load the eddy energy time series for the feedback run                                                                                                                               
    path_ee_fb = scratch_location+folder+run_fb+'output/energy.dat'
    ee_fb = np.loadtxt(path_ee_fb)
    ee_fb = ee_fb[:len(time),:]     #Keep only the relevant times                                                                                                                        

    feedback_ee = ee_fb[:,1]+ee_fb[:,2]

    loss[:,run_no+1] = control_ee - feedback_ee    #Absolute loss of energy


    for t in range(len(time)-1):

        power[t,run_no+1] = (loss[t+1,run_no+1]-loss[t,run_no+1])/(dt/s_to_day)


######################
# Make the damn plot #
######################

if(just_power==1):
    fig, ax2 = plt.subplots(figsize=(6,4))
    color = 'k'
    ax2.set_ylabel('Flow-to-wave energy flux (W/kg)', color=color)  # we already handled the x-label with ax1                                                                                                        
    ax2.tick_params(axis='y', labelcolor=color)
    #ax2.set_ylim(0, 10)                                                                                                                                                                                                  
    ax2.set_xlabel('Time (days)')
    ax2.grid(color=color, linestyle='-', linewidth=0.1)

     
    for run_no,run_fb in enumerate(run_list):
    
        ax2.plot(power[:,0],power[:,run_no+1],label='$u_0$ = '+storm_strength[run_no]+' cm/s',color=color,linestyle=linestyle_list[run_no],lw=1.5)

    ax2.legend(loc='best')

else:
    fig, ax1 = plt.subplots(figsize=(6,4))
    
    color = 'k'
    ax1.set_xlabel('Time (days)')
    ax1.set_ylabel('Total geostrophic energy loss', color=color)
    ax1.tick_params(axis='y', labelcolor=color)
    ax1.grid(color=color, linestyle='-', linewidth=0.1)
#    ax1.set_ylim(0, 2)
    
    for run_no,run_fb in enumerate(run_list):
        
        ax1.plot(time,loss[:,run_no+1],label='$u_0$ = '+storm_strength[run_no]+' cm/s',color=color,linestyle=linestyle_list[run_no],lw=1.5)
        
    plt.legend(loc='best',fontsize='small')
        
        
#Second axis: sfc EKE
    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis                                                                                                             
    
    color = 'b'
    ax2.set_ylabel('Flow-to-wave energy flux (W/kg)', color=color)  # we already handled the x-label with ax1    
    ax2.tick_params(axis='y', labelcolor=color)
#ax2.set_ylim(0, 10)
    
    for run_no,run_fb in enumerate(run_list):
        
        ax2.plot(power[:,0],power[:,run_no+1],color=color,linestyle=linestyle_list[run_no],lw=1.5)
        
fig.tight_layout()  # otherwise the right y-label is slightly clipped                                                                                 
plt.xlim(0,focus_time)


if (show==1):
    plt.show()                                                                                                                                                   
else:
    plt.savefig('plots/'+run+'/power.eps',bbox_inches='tight')















