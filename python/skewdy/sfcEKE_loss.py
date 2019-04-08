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

run = 'storm3/'    #Run to compare to                                                                                                                                                   
location = scratch_location+folder+run
 
run_list = ['storm3_fb/','storm3_fb_uw2/','storm3_fb_uw4/']
storm_strength = ['10','20','40']

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


#Now let's plot KE and WPE from the flow
path_ez = location+'output/ez.dat'
if os.path.isfile(path_ez):

    #Load the profiles time series                                                                                                                                                  
    ez = np.loadtxt(path_ez)

    #Extract Kinetic and potential energy of the flow
    fke = ez[:,1]  #Flow kinetic energy                                                                                                                        
    
    #Reshape so that the fields have dimensions wke[time,depth]                                                                                                                         
    fke = np.reshape(fke,(fke.shape[0]/n3,-1),order='C')

    #Keep only the focus region (in both depth and time)                                                                                                                                
    fke = fke[:max_time_ez,interest_depth]


    time = ts_ez_days*np.arange(0,fke.shape[0])

    loss = np.zeros( (len(time),len(run_list)+1) )
    loss[:,0] = time[:]

    
    #Import dats from the feedback runs
    for run_no,run_fb in enumerate(run_list):

        #Load the eddy energy time series for the feedback run                                                                                                                         
        path_ez_fb = scratch_location+folder+run_fb+'output/ez.dat'
        ez_fb = np.loadtxt(path_ez_fb)


        #Extract Kinetic and potential energy of the flow                                                                                                          
        fke_fb = ez_fb[:,1]  #Flow kinetic energy                                                                                                                             

        #Reshape so that the fields have dimensions wke[time,depth]                                                                                                         
        fke_fb = np.reshape(fke_fb,(fke_fb.shape[0]/n3,-1),order='C')

        #Keep only the focus region (in time, and keep only the surface)                                                                                                          
        fke_fb = fke_fb[:max_time_ez,interest_depth]

        
        loss[:,run_no+1] = 100*(fke - fke_fb)/fke  #Percentage of total eddy energy loss    


fig = plt.figure(figsize=(8,4))
ax = fig.add_subplot(1, 1, 1)
ax.grid(color='k', linestyle='-', linewidth=0.1)
for run_no,run_fb in enumerate(run_list):

    plt.plot(time,loss[:,run_no+1],label='$u_0$ = '+storm_strength[run_no]+' cm/s')

plt.legend(loc='best')
plt.xlabel('Time (days)')
plt.ylabel('Surface EKE loss (%)')
plt.xlim(0,focus_time)
#plt.show()                                                                                                                                                                              
plt.savefig('plots/'+run+'/sfcEKE_loss.eps',bbox_inches='tight')
