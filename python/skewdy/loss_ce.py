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
run    = 'wind_uw10/'    #Run to compare to                                                                                                                   
run_w  = 'wind_t30/'    #Run to compare to                               
run_wf = 'wind_uw10_t30/'    #Run to compare to  

location = scratch_location+folder+run    #Assumes the runs compared have the same parameters except feedback
location_w = scratch_location+folder+run_w    #Assumes the runs compared have the same parameters except feedback
location_wf = scratch_location+folder+run_wf    #Assumes the runs compared have the same parameters except feedback
 
focus_time =150  #Focus on the first $focus_time days  

show=1

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
path_ee  = location+'output/energy.dat'
path_ee_w  = location_w+'output/energy.dat'
path_ee_wf = location_wf+'output/energy.dat'

if os.path.isfile(path_ee) and os.path.isfile(path_ee_w) and os.path.isfile(path_ee_wf):

    #Load the profiles time series                                                                                                                                                  
    ee_f  = np.loadtxt(path_ee)
    ee_w  = np.loadtxt(path_ee_w)
    ee_wf = np.loadtxt(path_ee_wf)

    et_f  = ee_f[:max_time,1] + ee_f[:max_time,2]
    et_w  = ee_w[:max_time,1] + ee_w[:max_time,2]
    et_wf = ee_wf[:max_time,1] + ee_wf[:max_time,2]

    time = ts_e_days*np.arange(0,et_f.shape[0])

    loss_f  = 100*(et_w-et_f)/et_f[0]    #Loss percentage no wind (just feedback)
    loss_wf = 100*(et_w-et_wf)/et_f[0]    #Loss percentage with wind (no feedback)


#Load coupled energy
path_ce_f  = location+'output/ce.dat'
path_ce_wf = location_wf+'output/ce.dat'    

if os.path.isfile(path_ce_f) and os.path.isfile(path_ce_wf):

    #Load the profiles time series                                                                                                                                                  
    ce_f  = np.loadtxt(path_ce_f)
    ce_wf = np.loadtxt(path_ce_wf)

    ct_f  = 100*(ce_f[:max_time,1]  + ce_f[:max_time,2])/et_f[0]
    ct_wf = 100*(ce_wf[:max_time,1] + ce_wf[:max_time,2])/et_f[0]


######################
# Make the damn plot #
######################

fig, ax1 = plt.subplots(figsize=(6,4))

ax1.set_xlabel('Time (days)')
ax1.set_ylabel('Energy change (%)')
ax1.grid(color='k', linestyle='-', linewidth=0.1)
#ax1.set_ylim(0, 2)

ax1.plot(time,loss_f,label='No wind, EKE+EPE loss',color='b')
ax1.plot(time,ct_f,'--b',label='WPE+WCE')


ax1.plot(time,loss_wf,label=r'$\tau$ = 30 days, EKE+EPE loss',color='g')
ax1.plot(time,ct_wf,'--g',label=r'WPE+WCE')


plt.legend(loc='best',fontsize='small')


plt.xlim(0,focus_time)
if (show==1):
    plt.show()                                                                                                                                                   
else:
    plt.savefig('plots/'+run+'/loss_ce.eps',bbox_inches='tight')















