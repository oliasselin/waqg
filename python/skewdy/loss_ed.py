#!/usr/bin/env python                                                                                                                                                               
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker
from finds import find_resolution
from finds import find_scales
from finds import find_timestep


scratch_location = '/oasis/scratch/comet/oasselin/temp_project/'
folder = 'niskine/skewdy/'

u0='10'
run = 'wind/'    #Run to compare to                                                                                                                                                   
run_fb = 'wind_uw10/'    #Run to compare to                                                                                                                                            
location = scratch_location+folder+run    #Assumes the runs compared have the same parameters except feedback
location_fb = scratch_location+folder+run_fb    #Assumes the runs compared have the same parameters except feedback
 
focus_time =100  #Focus on the first $focus_time days  

show=0

#Read parameters from the source#
n1,n2,n3 = find_resolution(location)
Dx,Dz,L_scale,H_scale,U_scale,h_thermo = find_scales(location)
delt,freq_etot,freq_ez,freq_we,freq_wz = find_timestep(location)

T_scale = L_scale/U_scale
s_to_day=1./(3600*24)
ts_e_days = delt*freq_etot*T_scale*s_to_day
max_time = np.ceil(focus_time/ts_e_days)+1#number of time steps to reach focus_time days...


if not os.path.exists('plots/'+run):
    os.makedirs('plots/'+run)


#Now let's plot KE and WPE from the flow
path_ee = location+'output/energy.dat'
path_ee_fb = location_fb+'output/energy.dat'
if os.path.isfile(path_ee) and os.path.isfile(path_ee_fb):

    #Load the profiles time series                                                                                                                                                  
    ee = np.loadtxt(path_ee)
    ee_fb = np.loadtxt(path_ee_fb)

    et = ee[:max_time,1] + ee[:max_time,2]
    et_fb = ee_fb[:max_time,1] + ee_fb[:max_time,2]

    time = ts_e_days*np.arange(0,et.shape[0])

    loss = 100*(et-et_fb)/et[0]    #Loss percentage


######################
# Make the damn plot #
######################

fig, ax1 = plt.subplots(figsize=(6,4))

color = 'k'
ax1.set_xlabel('Time (days)')
ax1.set_ylabel('Total geostrophic energy loss (%)', color=color)
ax1.tick_params(axis='y', labelcolor=color)
ax1.grid(color=color, linestyle='-', linewidth=0.1)
#ax1.set_ylim(0, 2)

ax1.plot(time,loss,label='$u_0$ = '+u0+' cm/s',color=color)
#plt.legend(loc='best',fontsize='small')


plt.xlim(0,focus_time)
if (show==1):
    plt.show()                                                                                                                                                   
else:
    plt.savefig('plots/'+run+'/loss_ed.eps',bbox_inches='tight')















