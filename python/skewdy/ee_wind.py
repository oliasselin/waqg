#!/usr/bin/env python                                                                                                                                                               
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker
from finds import find_resolution
from finds import find_scales
from finds import find_timestep

tau_list=['','_t30','_t20','_t10']
label_list=['IVP','$ t = 30$ days','$ t = 20$ days','$ t = 10$ days']
line_list=['-b','-g','--g','.g']
u0='10'
run = 'wind_t30/'    #Control run: no feedback                     
run_root = 'wind_uw10'
focus_time =150  #Focus on the first $focus_time days  
show=0



scratch_location = '/oasis/scratch/comet/oasselin/temp_project/'
folder = 'niskine/skewdy/'
location      = scratch_location+folder+run    #Assumes the runs compared have the same parameters except feedback
location_root = scratch_location+folder+run_root    #Assumes the runs compared have the same parameters except feedback
 


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
path_ee    = location+'output/energy.dat'


et = np.zeros((max_time,len(tau_list)+1))
de = np.zeros((max_time,len(tau_list)+1))

if os.path.isfile(path_ee):


    #Load the geostrophic energy time series                        
    ee  = np.loadtxt(path_ee)
    et[:,0]  = ee[:max_time,1] + ee[:max_time,2]
    de[:,0]  = 100*(et[:,0]-et[0,0])/et[0,0]

    time = ts_e_days*np.arange(0,et.shape[0])

    for ntau,tau in enumerate(tau_list):

        eef = np.loadtxt(location_root+tau+'/output/energy.dat')
        et[:,ntau+1]  = eef[:max_time,1] + eef[:max_time,2]

        de[:,ntau+1] = 100*(et[:,ntau+1]-et[0,0])/et[0,0]

######################
# Make the damn plot #
######################

fig, ax2 = plt.subplots(figsize=(6,4))

ax2.set_xlabel('Time (days)')
ax2.set_ylabel('Total geostrophic energy change (%)')
ax2.grid(color='k', linestyle='-', linewidth=0.1)
#ax1.set_ylim(0, 2)

ax2.plot(time,de[:,0],'-k',label='No feedback')
for ntau,tau in enumerate(tau_list):
    ax2.plot(time,de[:,ntau+1],line_list[ntau],label=label_list[ntau])
    
xticks=np.arange(0,focus_time+1,30)
ax2.set_xticks(xticks)
ax2.legend(loc='best',fontsize='small')


plt.xlim(0,focus_time)
if (show==1):
    plt.show()                                                                                                                                                   
else:
    plt.savefig('plots/'+run+'/ee_wind.eps',bbox_inches='tight')















